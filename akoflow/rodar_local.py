"""
Executor local do `workflow.yaml` — para testar o DAG sem o AkôFlow.

Lê o mesmo arquivo que o AkôFlow lê e roda cada atividade como um processo,
respeitando `dependsOn` e um orçamento de CPU (soma dos `cpuLimit` em execução).
Diferença para o AkôFlow, declarada: aqui todas as atividades compartilham **um**
diretório de trabalho, em vez de um *workspace* por atividade com `rsync` das
predecessoras; o tempo medido não inclui o *overhead* de contêiner nem de `rsync`.

Uso::

    python akoflow/rodar_local.py --yaml akoflow/workflow.yaml --app-dir akoflow \\
        --workdir /tmp/G1/Zika_G1 --modo paralelo --tempos /tmp/G1/tempos.json
"""
import argparse
import json
import os
import subprocess
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor

import yaml

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _comando(atividade: dict, app_dir: str):
    cmd = atividade["command"]
    entrypoint = sys.executable if cmd["entrypoint"] == "python" else cmd["entrypoint"]
    args = [a.replace("/app/", app_dir.rstrip("/") + "/") if isinstance(a, str) else str(a)
            for a in cmd["arguments"]]
    return [entrypoint] + args


def _concluida(workdir: str, nome: str) -> bool:
    """Fragmento de proveniência com desfecho `concluido`: a atividade não precisa rodar de novo."""
    caminho = os.path.join(workdir, "out", "outputs", "prov", f"{nome}.json")
    try:
        with open(caminho, encoding="utf-8") as fh:
            return json.load(fh)["desfecho"]["status"] == "concluido"
    except (OSError, ValueError, KeyError):
        return False


def executar(yaml_path: str, app_dir: str, workdir: str, modo: str, jobs: int, tempos: str,
             retomar: bool = False) -> int:
    with open(yaml_path, encoding="utf-8") as fh:
        wf = yaml.safe_load(fh)
    atividades = {a["name"]: a for a in wf["spec"]["activities"]}
    cpu = {n: int(a.get("cpuLimit", 1)) for n, a in atividades.items()}
    orcamento = 1 if modo == "serial" else jobs
    os.makedirs(workdir, exist_ok=True)
    env = dict(os.environ, PYTHONPATH=REPO, MPLBACKEND="Agg")
    trava = threading.Condition()
    estado = {n: "ok" if retomar and _concluida(workdir, n) else "pendente" for n in atividades}
    emuso = [0]
    registro = {}
    t0 = time.time()

    def rodar(nome):
        with trava:
            trava.wait_for(lambda: emuso[0] + min(cpu[nome], orcamento) <= orcamento)
            emuso[0] += min(cpu[nome], orcamento)
        inicio = time.time()
        with open(os.path.join(workdir, f"{nome}.stdout.txt"), "w") as saida:
            codigo = subprocess.run(_comando(atividades[nome], app_dir), cwd=workdir, env=env,
                                    stdout=saida, stderr=subprocess.STDOUT).returncode
        fim = time.time()
        with trava:
            emuso[0] -= min(cpu[nome], orcamento)
            estado[nome] = "ok" if codigo == 0 else "falhou"
            registro[nome] = {"inicio_s": round(inicio - t0, 3), "fim_s": round(fim - t0, 3),
                              "duracao_s": round(fim - inicio, 3), "codigo": codigo, "cpu": cpu[nome]}
            trava.notify_all()

    with ThreadPoolExecutor(max_workers=max(1, jobs)) as pool:
        futuros = {}
        while True:
            with trava:
                prontos = [n for n, e in estado.items() if e == "pendente"
                           and all(estado[d] == "ok" for d in atividades[n].get("dependsOn", []))]
                bloqueadas = [n for n, e in estado.items() if e == "pendente"
                              and any(estado[d] in ("falhou", "bloqueada") for d in atividades[n].get("dependsOn", []))]
                for n in bloqueadas:
                    estado[n] = "bloqueada"
                for n in prontos:
                    estado[n] = "rodando"
                    futuros[n] = pool.submit(rodar, n)
                if all(e in ("ok", "falhou", "bloqueada") for e in estado.values()):
                    break
                trava.wait(timeout=0.2)
    total = time.time() - t0
    resumo = {"modo": modo, "jobs": orcamento, "makespan_s": round(total, 3),
              "trabalho_total_s": round(sum(r["duracao_s"] for r in registro.values()), 3),
              "atividades": registro, "bloqueadas": [n for n, e in estado.items() if e == "bloqueada"]}
    if tempos:
        with open(tempos, "w") as fh:
            json.dump(resumo, fh, indent=2)
    falhas = [n for n, e in estado.items() if e != "ok"]
    print(f"{modo}: makespan {total:.1f} s, {len(registro)} atividades, falhas/bloqueios: {falhas}")
    return 1 if falhas else 0


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--yaml", default=os.path.join(REPO, "akoflow", "workflow.yaml"))
    ap.add_argument("--app-dir", default=os.path.join(REPO, "akoflow"),
                    help="onde ficam config/ e data/ (o /app/ do YAML)")
    ap.add_argument("--workdir", required=True, help="diretório do projeto; o `out/` nasce nele")
    ap.add_argument("--modo", choices=("serial", "paralelo"), default="paralelo")
    ap.add_argument("--jobs", type=int, default=os.cpu_count(), help="orçamento de CPU no modo paralelo")
    ap.add_argument("--retomar", action="store_true",
                    help="pula as atividades cujo fragmento prov/ já diz `concluido`")
    ap.add_argument("--tempos", help="grava os tempos por atividade neste JSON")
    a = ap.parse_args()
    sys.exit(executar(a.yaml, os.path.abspath(a.app_dir), os.path.abspath(a.workdir), a.modo, a.jobs, a.tempos, a.retomar))
