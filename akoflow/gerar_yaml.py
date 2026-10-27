"""
Gera o `workflow.yaml` executável no AkôFlow local a partir do `akoflow/workflow.yaml` da G4b.

O YAML da G4b é o desenho (39 ativações, SUPOSIÇÕES marcadas). Este script o ajusta às
divergências de `docs/akoflow/07 §4` medidas no teste local (`docs/akoflow/08`):

* executável: `command.executable` (forma não legada) com a imagem local, no lugar de `spec.image`;
* `runtime`: o id do runtime do `environment.yaml` local (a dica da atividade não escolhe nada);
* conjunto de dados: configuração e diretório de entrada trocados (`--config`, `--entrada`);
* `expectedOutputs`: o fragmento de proveniência de toda atividade, mais as saídas de contrato;
* `dataDependencies` ⊂ `dependsOn`: conferido, e `sizeBytes` do alinhamento reescalado.

Uso::

    python akoflow/gerar_yaml.py --nome ptm-zika21 --config /app/config/zika21.json \\
        --entrada /app/data/Zika479_Test_large --saida /tmp/zika21.yaml
"""
import argparse
import sys

import yaml

ZIKA_CONFIG, ZIKA_ENTRADA = "/app/config/zika21.json", "/app/data/Zika479_Test_large"
ZIKA_GB = ZIKA_ENTRADA + "/raw_data_sequences.gb"


def _trocar(args, antigo, novo):
    return [a.replace(antigo, novo) if isinstance(a, str) else a for a in args]


def gerar(base: dict, nome: str, config: str, entrada: str, imagem: str, runtime: str,
          namespace: str, tamanho_aln: int = None, sem_gb: bool = False,
          serializar_saidas: bool = False, empacotar: bool = False) -> dict:
    wf = yaml.safe_load(yaml.safe_dump(base))
    spec = wf["spec"]
    spec.pop("image", None)  # forma legada: vale `command.executable` por atividade
    spec["namespace"] = namespace
    wf["name"] = nome
    nomes = {a["name"] for a in spec["activities"]}
    for dep in spec.get("dataDependencies", []):
        produtor = next(a for a in spec["activities"] if a["name"] == dep["consumerActivity"])
        if dep["producerActivity"] not in produtor.get("dependsOn", []):
            raise SystemExit(f"dataDependency sem dependsOn: {dep['producerActivity']} -> {dep['consumerActivity']}")
        if tamanho_aln and dep["logicalName"].endswith(".aln"):
            dep["sizeBytes"] = tamanho_aln
    for a in spec["activities"]:
        cmd = a["command"]
        cmd["arguments"] = _trocar(_trocar(cmd["arguments"], ZIKA_CONFIG, config), ZIKA_ENTRADA, entrada)
        cmd["executable"] = {"source": {"type": "oci", "reference": imagem}, "delivery": {"strategy": "auto"}}
        a["runtime"] = runtime
        # Todo fragmento de proveniência existe, mesmo quando a ativação é no-op (D28): é o contrato mínimo.
        # Exceto `consolidar-manifesto`: ele LÊ os fragmentos e escreve o manifesto, não grava fragmento
        # próprio (medido: o AkôFlow cancela a atividade se um expectedOutput não for observado).
        saidas = [] if a["name"] == "consolidar-manifesto" else [f"out/outputs/prov/{a['name']}.json"]
        if a["name"] == "preparar-entrada":
            saidas += ["out/outputs/config_backup.json", "out/entrada/acessos.txt"]
        if a["name"] == "consolidar-manifesto":
            saidas += ["out/outputs/manifest.json", "out/outputs/output_log.txt"]
        a["expectedOutputs"] = saidas
    if sem_gb:  # sem `.gb` congelado na imagem (VARV-6): `baixar-metadados` baixa do NCBI (modo online)
        for a in spec["activities"]:
            args = a["command"]["arguments"]
            if "--gb-congelado" in args:
                i = args.index("--gb-congelado")
                del args[i:i + 2]
    if serializar_saidas:  # contorno da corrida do daemon (08): dois consumidores de `minerar-fpmax` ao mesmo tempo
        graficos = next(a for a in spec["activities"] if a["name"] == "graficos-suporte")
        if "gerar-cql" not in graficos["dependsOn"]:
            graficos["dependsOn"] = list(graficos["dependsOn"]) + ["gerar-cql"]
    if empacotar:  # o daemon poda as entradas herdadas do workspace final (08): só o que a última atividade GRAVA sobrevive
        spec["activities"].append({
            "name": "empacotar-saida", "runtime": runtime, "cpuLimit": "1", "memoryLimit": "1Gi",
            "dependsOn": ["consolidar-manifesto"],
            "command": {"entrypoint": "sh",
                        "arguments": ["-c", "mkdir -p resultados && tar -cf resultados/out.tar out"],
                        "executable": {"source": {"type": "oci", "reference": imagem}, "delivery": {"strategy": "auto"}}},
            "expectedOutputs": ["resultados/out.tar"]})
    assert len(nomes) + int(empacotar) == len(spec["activities"])
    return wf


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--base", default="akoflow/workflow.yaml")
    ap.add_argument("--nome", required=True)
    ap.add_argument("--config", default=ZIKA_CONFIG)
    ap.add_argument("--entrada", default=ZIKA_ENTRADA)
    ap.add_argument("--imagem", default="phylotreeminer-akoflow:local")
    ap.add_argument("--runtime", default="ptm-local-runtime")
    ap.add_argument("--namespace", default="ptm")
    ap.add_argument("--tamanho-aln", type=int, help="bytes de um alinhamento (sizeBytes das dataDependencies)")
    ap.add_argument("--sem-gb", action="store_true", help="remove --gb-congelado (baixa do NCBI)")
    ap.add_argument("--serializar-saidas", action="store_true",
                    help="graficos-suporte depende também de gerar-cql (evita a corrida de rsync do daemon)")
    ap.add_argument("--empacotar", action="store_true",
                    help="acrescenta `empacotar-saida`: tar do out/ herdado, que o daemon de outro modo poda")
    ap.add_argument("--saida", required=True)
    a = ap.parse_args()
    with open(a.base, encoding="utf-8") as fh:
        base = yaml.safe_load(fh)
    wf = gerar(base, a.nome, a.config, a.entrada, a.imagem, a.runtime, a.namespace, a.tamanho_aln, a.sem_gb, a.serializar_saidas, a.empacotar)
    with open(a.saida, "w", encoding="utf-8") as fh:
        yaml.safe_dump(wf, fh, sort_keys=False, allow_unicode=True)
    print(f"{a.saida}: {len(wf['spec']['activities'])} atividades", file=sys.stderr)
