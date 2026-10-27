"""
Registra um workflow no AkôFlow local, planeja (HEFT), executa e coleta o que o servidor gravou.

Fluxo da API (docs/akoflow/07 §2.7-2.9, 2.18, medido em 08): ambiente -> escopo -> topologia ->
definição do workflow -> sessão de planejamento HEFT -> candidato selecionado -> POST
/execution-runs/ com o *bundle* (run + plano + workflow + recursos + escopo + topologia + runtimes).
O corpo só com `run` falha com `FOREIGN KEY constraint failed` (08 §problemas).

Uso::

    export AKOFLOW_API_URL=http://127.0.0.1:8080/akoflow-api AKOFLOW_API_TOKEN=...
    python akoflow/rodar_akoflow.py --yaml zika21.yaml --ambiente akoflow/environment.local.yaml \\
        --id zika21-r1 --coleta /tmp/coleta/zika21-r1
"""
import argparse
import json
import os
import sys
import time
import urllib.error
import urllib.request

import yaml

API = os.environ.get("AKOFLOW_API_URL", "http://127.0.0.1:8080/akoflow-api").rstrip("/")
TOKEN = os.environ["AKOFLOW_API_TOKEN"]


def chamar(metodo, caminho, corpo=None, tipo="application/json", tolerar=()):
    dados = None
    if corpo is not None:
        dados = corpo if isinstance(corpo, bytes) else json.dumps(corpo).encode()
    req = urllib.request.Request(f"{API}/{caminho.lstrip('/')}", data=dados, method=metodo,
                                 headers={"Authorization": f"Bearer {TOKEN}", "Content-Type": tipo})
    try:
        with urllib.request.urlopen(req, timeout=120) as r:
            texto = r.read().decode()
    except urllib.error.HTTPError as e:
        texto = e.read().decode()
        if e.code in tolerar:
            return {"_http": e.code, "_corpo": texto}
        raise SystemExit(f"{metodo} {caminho} -> HTTP {e.code}: {texto[:600]}")
    if not texto.strip():
        return {}
    try:
        return json.loads(texto)
    except ValueError:  # alguns endpoints de proveniência devolvem texto
        return {"_texto": texto}


def registrar_ambiente(arq):
    with open(arq, "rb") as fh:
        bruto = fh.read()
    amb = yaml.safe_load(bruto)
    ver = amb["version"]["id"]
    r = chamar("POST", "environments/", bruto, "application/yaml", tolerar=(409, 422, 500))
    escopo = {"id": f"{amb['environment']['id']}-scope", "name": "escopo local", "environmentVersionIds": [ver]}
    chamar("POST", "execution-scopes/", escopo, tolerar=(409, 422, 500))
    topo = {"id": f"{amb['environment']['id']}-topology", "name": "topologia local", "version": 1,
            "executionScopeId": escopo["id"], "links": []}
    chamar("POST", "network-topologies/", topo, tolerar=(409, 422, 500))
    return amb, escopo["id"], topo["id"]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--yaml")
    ap.add_argument("--ambiente")
    ap.add_argument("--id", required=True, help="id do run (único)")
    ap.add_argument("--coleta", required=True)
    ap.add_argument("--coletar", action="store_true", help="não submete: só espera o run `--id` já criado e coleta")
    ap.add_argument("--espera-max", type=int, default=7200)
    a = ap.parse_args()
    os.makedirs(a.coleta, exist_ok=True)
    plano = {}
    if not a.coletar:
        amb, escopo, topo = registrar_ambiente(a.ambiente)
        env_id = amb["environment"]["id"]
        with open(a.yaml, "rb") as fh:
            yaml_bruto = fh.read()
        wf = chamar("POST", "workflow-definitions/", yaml_bruto, "application/yaml", tolerar=(422, 500))
        if "_http" in wf:  # `UNIQUE constraint failed: external_id`: a definição já existe; reaproveita (reexecução)
            wf = chamar("GET", f"workflow-definitions/{yaml.safe_load(yaml_bruto)['name']}/")
        wv = wf["version"]["id"]
        ps = f"{a.id}-plan"
        chamar("POST", "planning-sessions/", {"id": ps, "workflowVersionId": wv, "executionScopeId": escopo,
               "networkTopologyId": topo, "algorithms": [{"id": "heft", "configuration": {}}],
               "deadlineSeconds": 0, "budget": 0})
        for _ in range(120):
            s = chamar("GET", f"planning-sessions/{ps}/")
            if s["session"]["status"] in ("completed", "failed"):
                break
            time.sleep(1)
        cands = chamar("GET", f"planning-sessions/{ps}/candidates/")
        if not cands:
            raise SystemExit(f"sem candidato de plano: {s}")
        plano = chamar("POST", f"planning-sessions/{ps}/candidates/{cands[0]['id']}/select/")
        bundle = {"run": {"id": a.id, "schedulePlanId": plano["id"], "mode": "real", "status": "created"},
                  "plan": plano, "workflow": wf["version"],
                  "resources": chamar("GET", f"environments/{env_id}/")["resources"],
                  "executionScope": chamar("GET", f"execution-scopes/{escopo}/"),
                  "networkTopology": chamar("GET", f"network-topologies/{topo}/"),
                  "runtimes": chamar("GET", f"environments/{env_id}/")["runtimes"],
                  "runtimeBindings": chamar("GET", f"environments/{env_id}/").get("resourceRuntimeBindings", [])}
        chamar("POST", "execution-runs/", bundle)
    t0 = time.time()
    while time.time() - t0 < a.espera_max:
        r = chamar("GET", f"execution-runs/{a.id}/", tolerar=(404,))
        if "run" in r and r["run"]["status"] in ("completed", "failed", "cancelled"):
            break
        time.sleep(2)
    fim = time.time() - t0
    r = chamar("GET", f"execution-runs/{a.id}/")
    json.dump(r, open(os.path.join(a.coleta, "run.json"), "w"), indent=1)
    json.dump(plano, open(os.path.join(a.coleta, "plano.json"), "w"), indent=1)
    for nome, caminho in (("linhagem", f"provenance/lineage/runs/{a.id}/?direction=both&depth=3&maxNodes=300"),
                          ("entidades", "provenance/entities/"), ("auditoria", "audit-events/?limit=200")):
        try:
            json.dump(chamar("GET", caminho), open(os.path.join(a.coleta, f"{nome}.json"), "w"), indent=1)
        except (SystemExit, OSError) as e:
            print(f"[coleta {nome}] {e}", file=sys.stderr)
    print(f"{a.id}: {r['run']['status']} em {fim:.1f} s (relógio de parede, inclui o polling)")
    if r["run"].get("failureReason"):
        print("falha:", r["run"]["failureReason"][:500])
    return 0 if r["run"]["status"] == "completed" else 1


if __name__ == "__main__":
    sys.exit(main())
