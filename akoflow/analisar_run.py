"""
Resume um run do AkôFlow a partir do SQL de proveniência do próprio daemon (`POST /provenance/sql/`).

Saída: tabela atividade x início relativo x duração x estado, e as métricas do run: makespan,
trabalho total (soma das durações), caminho crítico medido (DAG de `activity_dependencies`),
paralelismo máximo e médio observados, e o speedup contra um tempo de monolito dado.

Uso::

    export AKOFLOW_API_URL=http://127.0.0.1:8080/akoflow-api AKOFLOW_API_TOKEN=...
    python akoflow/analisar_run.py zika21-r5 --mono 672 [--json saida.json]
"""
import argparse
import json
import os
import urllib.request

API = os.environ.get("AKOFLOW_API_URL", "http://127.0.0.1:8080/akoflow-api").rstrip("/")
TOKEN = os.environ["AKOFLOW_API_TOKEN"]


def sql(consulta: str, **parametros):
    corpo = json.dumps({"sql": consulta, "parameters": parametros, "page": 1, "pageSize": 200}).encode()
    req = urllib.request.Request(f"{API}/provenance/sql/", data=corpo, method="POST",
                                 headers={"Authorization": f"Bearer {TOKEN}", "Content-Type": "application/json"})
    with urllib.request.urlopen(req, timeout=60) as r:
        return json.loads(r.read())["items"]


def analisar(run: str, mono: float = None) -> dict:
    ex = sql("SELECT started_at, finished_at, status, failure_reason FROM execution_runs WHERE id = :r", r=run)[0]
    tarefas = sql("SELECT activity_id, status, started_at, finished_at, runtime_seconds, queue_seconds, "
                  "overhead_seconds FROM task_executions WHERE execution_run_id = :r ORDER BY started_at", r=run)
    deps = sql("SELECT d.activity_id a, d.depends_on_activity_id b FROM activity_dependencies d "
               "WHERE d.activity_id IN (SELECT activity_id FROM task_executions WHERE execution_run_id = :r)", r=run)
    t = {x["activity_id"]: x for x in tarefas if x["started_at"]}
    t0 = min(x["started_at"] for x in t.values())
    t1 = max(x["finished_at"] for x in t.values())
    makespan = t1 - t0
    # caminho crítico medido: maior soma de durações ao longo do DAG
    pred = {}
    for d in deps:
        pred.setdefault(d["a"], []).append(d["b"])
    memo = {}

    def cp(a):
        if a not in memo:
            ps = [p for p in pred.get(a, []) if p in t]
            melhor = max(ps, key=cp, default=None)
            memo[a] = (t[a]["runtime_seconds"] + (cp(melhor)[0] if melhor else 0), melhor)
        return memo[a]

    fim = max(t, key=lambda a: cp(a)[0])
    caminho, a = [], fim
    while a:
        caminho.append(a)
        a = cp(a)[1]
    eventos = sorted([(x["started_at"], 1) for x in t.values()] + [(x["finished_at"], -1) for x in t.values()])
    atual = maximo = 0
    for _, delta in eventos:
        atual += delta
        maximo = max(maximo, atual)
    trabalho = sum(x["runtime_seconds"] for x in t.values())
    r = {"run": run, "status": ex["status"], "falha": ex["failure_reason"], "atividades": len(t),
         "makespan_s": makespan, "trabalho_total_s": trabalho,
         "caminho_critico_s": cp(fim)[0], "caminho": [c.split("-", 2)[-1] for c in reversed(caminho)],
         "paralelismo_max": maximo, "paralelismo_medio": trabalho / makespan,
         "tarefas": [{"atividade": k, "inicio_s": v["started_at"] - t0, "duracao_s": v["runtime_seconds"],
                      "fila_s": v["queue_seconds"], "estado": v["status"]} for k, v in t.items()]}
    if mono:
        r["speedup_vs_mono"] = mono / makespan
    return r


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run")
    ap.add_argument("--mono", type=float, help="makespan do monolito (s)")
    ap.add_argument("--json")
    a = ap.parse_args()
    r = analisar(a.run, a.mono)
    for k in ("run", "status", "atividades", "makespan_s", "trabalho_total_s", "caminho_critico_s",
              "paralelismo_max", "paralelismo_medio", "speedup_vs_mono"):
        print(f"{k}: {r.get(k)}")
    print("caminho:", " > ".join(r["caminho"]))
    if r["falha"]:
        print("falha:", r["falha"][:300])
    for x in r["tarefas"]:
        print(f"{x['atividade'].split('-', 2)[-1]:50s} +{x['inicio_s']:8.1f}  {x['duracao_s']:8.2f}  {x['estado']}")
    if a.json:
        json.dump(r, open(a.json, "w"), indent=1)
