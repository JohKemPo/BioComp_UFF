"""Modelo do DAG granular (G4a): ativações, largura, caminho crítico e esboço YAML.

Uso: python dag.py zika|varv [--yaml] [--mermaid]
Durações: segundos, medidas nos logs/manifestos citados em 02-desenho-granular.md.
"""
import sys, json
import networkx as nx
from networkx.algorithms import bipartite

CFG = {
    "zika": dict(stem="dataset_final", aligners=["mafft", "mafft_iterative"],
                 basicos=["nj_distance", "upgma_distance", "nj_parsimony", "upgma_parsimony"],
                 avancados=["fasttree", "iqtree", "raxml"], align_threads=4, raxml_threads=4,
                 falham={("mafft", m) for m in ["nj_distance", "upgma_distance", "nj_parsimony", "upgma_parsimony"]},
                 T_mono=931.566),
    "varv": dict(stem="dataset_final", aligners=["mafft", "mafft_iterative"],
                 basicos=["nj_distance", "upgma_distance"],
                 avancados=["fasttree", "iqtree", "raxml"], align_threads=16, raxml_threads=8,
                 falham=set(), T_mono=9039.430),
}

DUR = {
    "zika": {
        "preparar-entrada": 0.003, "baixar-metadados": 4.508,
        "alinhar-mafft": 6.086, "alinhar-mafft-iterative": 10.06,
        "modelo-mafft": 1.256, "modelo-mafft-iterative": 1.152,
        "arvore-mafft-fasttree": 2.77 + 0.28, "arvore-mafft-iqtree": 6.216 + 0.43, "arvore-mafft-raxml": 271.348 + 0.26,
        "arvore-mafft-iterative-nj-distance": 0.11 + 0.22, "arvore-mafft-iterative-upgma-distance": 0.111 + 0.22,
        "arvore-mafft-iterative-nj-parsimony": 246.353 + 0.23, "arvore-mafft-iterative-upgma-parsimony": 78.789 + 0.22,
        "arvore-mafft-iterative-fasttree": 2.97 + 0.28, "arvore-mafft-iterative-iqtree": 5.141 + 0.32,
        "arvore-mafft-iterative-raxml": 285.481 + 0.28,
        "_falha": 0.01, "_subarvore": 0.42, "heatmap-rf": 1.6, "minerar-fpmax": 0.726,
        "gerar-cql": 0.0, "graficos-suporte": 0.0, "consolidar-manifesto": 0.0,
    },
    "varv": {
        "preparar-entrada": 0.085, "baixar-metadados": 35.047,
        "alinhar-mafft": 602.7, "alinhar-mafft-iterative": 2674.2,
        "modelo-mafft": 9.947, "modelo-mafft-iterative": 9.787,
        "arvore-mafft-nj-distance": 25.7, "arvore-mafft-upgma-distance": 25.691,
        "arvore-mafft-fasttree": 233.79, "arvore-mafft-iqtree": 46.673, "arvore-mafft-raxml": 2129.041,
        "arvore-mafft-iterative-nj-distance": 25.7, "arvore-mafft-iterative-upgma-distance": 25.651,
        "arvore-mafft-iterative-fasttree": 236.97, "arvore-mafft-iterative-iqtree": 46.747,
        "arvore-mafft-iterative-raxml": 2039.65,
        "_falha": 0.0, "_subarvore": 334.457, "heatmap-rf": 2.4, "minerar-fpmax": 20.234,
        "gerar-cql": 0.0, "graficos-suporte": 0.0, "consolidar-manifesto": 0.0,
    },
}

def h(s):
    return s.replace("_", "-")

def build(nome):
    c = CFG[nome]
    G = nx.DiGraph()
    def add(n, **kw):
        G.add_node(n, **kw)
    add("preparar-entrada", cpu=1)
    add("baixar-metadados", cpu=1); G.add_edge("preparar-entrada", "baixar-metadados")
    arvores = []
    for a in c["aligners"]:
        al = f"alinhar-{h(a)}"; add(al, cpu=c["align_threads"]); G.add_edge("preparar-entrada", al)
        mo = f"modelo-{h(a)}"; add(mo, cpu=1); G.add_edge(al, mo)
        for m in c["basicos"] + c["avancados"]:
            ar = f"arvore-{h(a)}-{h(m)}"
            add(ar, cpu=(c["raxml_threads"] if m == "raxml" else 1), aligner=a, metodo=m,
                falha=(a, m) in c["falham"])
            G.add_edge(al, ar)
            if m in ("iqtree", "raxml", "mrbayes"):
                G.add_edge(mo, ar)
            sb = f"subarvores-{h(a)}-{h(m)}"; add(sb, cpu=1)
            G.add_edge(ar, sb); G.add_edge("baixar-metadados", sb)
            arvores.append(ar)
    add("heatmap-rf", cpu=1)
    for ar in arvores: G.add_edge(ar, "heatmap-rf")
    add("minerar-fpmax", cpu=1)
    for n in list(G.nodes):
        if n.startswith("subarvores-"): G.add_edge(n, "minerar-fpmax")
    add("gerar-cql", cpu=1); G.add_edge("minerar-fpmax", "gerar-cql")
    add("graficos-suporte", cpu=1); G.add_edge("minerar-fpmax", "graficos-suporte")
    add("consolidar-manifesto", cpu=1)
    for n in ("heatmap-rf", "gerar-cql", "graficos-suporte"): G.add_edge(n, "consolidar-manifesto")
    assert nx.is_directed_acyclic_graph(G)
    d = DUR[nome]
    for n in G.nodes:
        if n in d: G.nodes[n]["dur"] = d[n]
        elif n.startswith("arvore-"): G.nodes[n]["dur"] = d["_falha"] if G.nodes[n].get("falha") else None
        elif n.startswith("subarvores-"):
            ar = "arvore-" + n[len("subarvores-"):]
            G.nodes[n]["dur"] = 0.0 if G.nodes[ar].get("falha") else d["_subarvore"]
    faltam = [n for n in G.nodes if G.nodes[n].get("dur") is None]
    assert not faltam, faltam
    return G

def largura(G):
    # Dilworth: maior anticadeia = n - emparelhamento máximo no fecho transitivo.
    T = nx.transitive_closure_dag(G)
    B = nx.Graph()
    L = [("L", n) for n in G]; R = [("R", n) for n in G]
    B.add_nodes_from(L, bipartite=0); B.add_nodes_from(R, bipartite=1)
    B.add_edges_from((("L", u), ("R", v)) for u, v in T.edges)
    M = bipartite.hopcroft_karp_matching(B, top_nodes=L)
    emparelhados = sum(1 for k in M if k[0] == "L")
    return G.number_of_nodes() - emparelhados

def caminho_critico(G):
    fim = {}; pred = {}
    for n in nx.topological_sort(G):
        ini = max((fim[p] for p in G.predecessors(n)), default=0.0)
        pred[n] = max(G.predecessors(n), key=lambda p: fim[p], default=None)
        fim[n] = ini + G.nodes[n]["dur"]
    ult = max(fim, key=fim.get)
    cam = []
    while ult: cam.append(ult); ult = pred[ult]
    return fim, list(reversed(cam))

def pico_asap(G, fim):
    ev = []
    for n in G:
        f = fim[n]; i = f - G.nodes[n]["dur"]
        if G.nodes[n]["dur"] > 0: ev += [(i, 1, G.nodes[n]["cpu"]), (f, -1, -G.nodes[n]["cpu"])]
    ev.sort(key=lambda e: (e[0], e[1]))
    a = c = pa = pc = 0
    for _, da, dc in ev:
        a += da; c += dc; pa = max(pa, a); pc = max(pc, c)
    return pa, pc

def resumo():
    nome = sys.argv[1]
    G = build(nome)
    fim, cam = caminho_critico(G)
    cp = max(fim.values())
    W = sum(G.nodes[n]["dur"] for n in G)
    pa, pc = pico_asap(G, fim)
    print(json.dumps({
        "ativacoes": G.number_of_nodes(), "arestas": G.number_of_edges(),
        "por_tipo": {k: sum(1 for n in G if n.startswith(k)) for k in
                     ["preparar", "baixar", "alinhar", "modelo", "arvore", "subarvores", "heatmap", "minerar", "gerar", "graficos", "consolidar"]},
        "largura_max_anticadeia": largura(G),
        "caminho_critico_s": round(cp, 1), "caminho": [(n, round(G.nodes[n]["dur"], 3)) for n in cam],
        "trabalho_total_s": round(W, 1), "T_mono_s": CFG[nome]["T_mono"],
        "speedup_teorico_vs_mono": round(CFG[nome]["T_mono"] / cp, 2),
        "paralelismo_medio_W_sobre_CP": round(W / cp, 2),
        "pico_ativacoes_simultaneas_ASAP": pa, "pico_cpu_ASAP": pc,
    }, indent=1, ensure_ascii=False))

if __name__ == "__main__" and len(sys.argv) == 2:
    resumo()

# ---------------------------------------------------------------- YAML
MEM = {"preparar-entrada": "512Mi", "baixar-metadados": "512Mi", "alinhar": "2Gi", "modelo": "1Gi",
       "arvore": "2Gi", "subarvores": "1Gi", "heatmap-rf": "1Gi", "minerar-fpmax": "2Gi",
       "gerar-cql": "1Gi", "graficos-suporte": "1Gi", "consolidar-manifesto": "1Gi"}
SIZES = {"zika": {"aln": 220915}, "varv": {"aln": 11120642}}

def mem(n):
    for k, v in MEM.items():
        if n == k or n.startswith(k + "-") or n.startswith(k):
            return v
    raise KeyError(n)

def comando(G, n, nome):
    c = CFG[nome]; cfgp = "out/outputs/config_backup.json"
    base = ["-m", "workflow.akoflow"]
    if n == "preparar-entrada":
        return base + ["preparar", "--config", "/app/config/zika21.json", "--entrada", "/app/data/Zika479_Test_large", "--out", "out"]
    if n == "baixar-metadados":
        return base + ["baixar-metadados", "--config", cfgp, "--acessos", "out/entrada/acessos.txt",
                       "--gb-congelado", "/app/data/Zika479_Test_large/raw_data_sequences.gb", "--out", "out"]
    if n.startswith("alinhar-"):
        a = n[len("alinhar-"):].replace("-", "_")
        return base + ["alinhar", "--config", cfgp, "--alinhador", a, "--out", "out"]
    if n.startswith("modelo-"):
        a = n[len("modelo-"):].replace("-", "_")
        return base + ["selecionar-modelo", "--config", cfgp, "--alinhador", a, "--out", "out"]
    if n.startswith("arvore-"):
        nd = G.nodes[n]
        return base + ["inferir", "--config", cfgp, "--alinhador", nd["aligner"], "--metodo", nd["metodo"], "--out", "out"]
    if n.startswith("subarvores-"):
        ar = G.nodes["arvore-" + n[len("subarvores-"):]]
        arv = f"out/Trees/tree_{c['stem']}_{ar['aligner']}_{ar['metodo']}.nexus"
        return base + ["subarvores", "--config", cfgp, "--arvore", arv, "--out", "out"]
    if n == "heatmap-rf":
        return base + ["heatmap-rf", "--config", cfgp, "--stem", c["stem"], "--out", "out"]
    if n == "minerar-fpmax":
        return base + ["minerar", "--config", cfgp, "--out", "out"]
    if n == "gerar-cql":
        return base + ["gerar-cql", "--config", cfgp, "--out", "out"]
    if n == "graficos-suporte":
        return base + ["graficos-suporte", "--config", cfgp, "--out", "out"]
    if n == "consolidar-manifesto":
        return base + ["consolidar", "--config", cfgp, "--out", "out"]
    raise KeyError(n)

def yaml_esboco(nome):
    G = build(nome)
    L = []
    L.append("name: phylotreeminer-zika21-granular")
    L.append("spec:")
    L.append("  activities:")
    for n in nx.lexicographical_topological_sort(G, key=lambda x: x):
        nd = G.nodes[n]
        L.append("    - command:")
        L.append("        arguments:")
        for a in comando(G, n, nome):
            L.append(f"          - {a}")
        L.append("        entrypoint: python")
        L.append(f"      cpuLimit: \"{nd['cpu']}\"")
        deps = sorted(G.predecessors(n))
        if deps:
            L.append("      dependsOn:")
            for d in deps:
                L.append(f"        - {d}")
        L.append(f"      memoryLimit: {mem(n)}  # SUPOSIÇÃO (valor não medido)")
        L.append(f"      name: {n}")
        L.append("      runtime: phylotreeminer-local  # SUPOSIÇÃO (id definido no environment.yaml da equipe)")
    L.append("  dataDependencies:  # campo do exemplo montage-050 (mesmo commit), ausente no infrastructure-selection")
    for n in sorted(G.nodes):
        if n.startswith("arvore-") or n.startswith("modelo-"):
            for p in sorted(G.predecessors(n)):
                if p.startswith("alinhar-"):
                    a = p[len("alinhar-"):].replace("-", "_")
                    L.append(f"    - consumerActivity: {n}")
                    L.append(f"      logicalName: tmp/{CFG[nome]['stem']}_{a}.aln")
                    L.append(f"      producerActivity: {p}")
                    L.append(f"      sizeBytes: {SIZES[nome]['aln']}  # medido: tamanho do .aln no projeto de referência")
    L.append("  image: ghcr.io/SUPOSICAO/phylotreeminer-workflow:29833ff  # SUPOSIÇÃO: registro e tag a definir com a equipe")
    L.append("  namespace: phylotreeminer  # SUPOSIÇÃO: o exemplo usa 'showcase'")
    return "\n".join(L) + "\n"

# ------------------------------------------------------------- Mermaid
def mid(n):
    return n.replace("-", "_")

def mermaid_expandido(nome):
    G = build(nome); c = CFG[nome]; s = c["stem"]
    L = ["flowchart LR"]
    L.append('  F_in[/"data/Zika479_Test_large/dataset_final.fasta<br/>55 820 B"/]')
    L.append('  F_cfg[/"config zika21.json"/]')
    L.append('  F_gbref[/"raw_data_sequences.gb congelado<br/>(do projeto de referência)"/]')
    def act(n, extra=""):
        rot = n + extra
        return f'  {mid(n)}["{rot}"]'
    L.append(act("preparar-entrada"))
    L.append("  F_in --> preparar_entrada"); L.append("  F_cfg --> preparar_entrada")
    L.append('  F_fa[/"out/entrada/dataset_final.fasta<br/>+ acessos.txt + config_backup.json"/]')
    L.append("  preparar_entrada --> F_fa")
    L.append(act("baixar-metadados")); L.append("  F_fa --> baixar_metadados"); L.append("  F_gbref --> baixar_metadados")
    L.append('  F_gb[/"out/outputs/raw_data_sequences.gb"/]'); L.append("  baixar_metadados --> F_gb")
    for a in c["aligners"]:
        ah = h(a)
        L.append(f'  subgraph S_{mid(ah)}["braço {a}"]')
        L.append("    direction LR")
        L.append(f'    {mid("alinhar-"+ah)}["alinhar-{ah}<br/>cpu {c["align_threads"]}"]')
        L.append(f'    F_aln_{mid(ah)}[/"tmp/{s}_{a}.aln<br/>Align/{s}_{a}.aln"/]')
        L.append(f'    {mid("modelo-"+ah)}["modelo-{ah}<br/>ModelFinder BIC"]')
        L.append(f'    F_mod_{mid(ah)}[/"tmp/modelfinder_{s}_{a}/<br/>selecao.json"/]')
        for m in c["basicos"] + c["avancados"]:
            ar = f"arvore-{ah}-{h(m)}"; sb = f"subarvores-{ah}-{h(m)}"
            falha = G.nodes[ar].get("falha")
            rot = ar + ("<br/>tentado_e_falhou (D28)" if falha else "")
            L.append(f'    {mid(ar)}["{rot}"]')
            if falha:
                L.append(f'    F_t_{mid(ar)}[/"sem árvore<br/>prov/{ar}.json"/]')
            else:
                L.append(f'    F_t_{mid(ar)}[/"Trees/tree_{s}_{a}_{m}.nexus<br/>Plots/*.png"/]')
            L.append(f'    {mid(sb)}["{sb}"]')
            L.append(f'    F_s_{mid(sb)}[/"Subtrees/*_InnerN.nexus<br/>outputs/subarvores/*.json"/]')
        L.append("  end")
        L.append(f"  F_fa --> {mid('alinhar-'+ah)}")
        L.append(f"  {mid('alinhar-'+ah)} --> F_aln_{mid(ah)}")
        L.append(f"  F_aln_{mid(ah)} --> {mid('modelo-'+ah)}")
        L.append(f"  {mid('modelo-'+ah)} --> F_mod_{mid(ah)}")
        for m in c["basicos"] + c["avancados"]:
            ar = f"arvore-{ah}-{h(m)}"; sb = f"subarvores-{ah}-{h(m)}"
            L.append(f"  F_aln_{mid(ah)} --> {mid(ar)}")
            if m in ("iqtree", "raxml", "mrbayes"):
                L.append(f"  F_mod_{mid(ah)} --> {mid(ar)}")
            L.append(f"  {mid(ar)} --> F_t_{mid(ar)}")
            L.append(f"  F_t_{mid(ar)} --> {mid(sb)}")
            L.append(f"  F_gb --> {mid(sb)}")
            L.append(f"  {mid(sb)} --> F_s_{mid(sb)}")
            L.append(f"  F_s_{mid(sb)} --> minerar_fpmax")
            L.append(f"  F_t_{mid(ar)} --> heatmap_rf")
    L.append('  heatmap_rf["heatmap-rf"]'); L.append('  F_hm[/"Plots/tree_dataset_final_heatmap_of_distances.png"/]')
    L.append("  heatmap_rf --> F_hm")
    L.append('  minerar_fpmax["minerar-fpmax<br/>FPMax OFST, suporte auto"]')
    L.append('  F_meta[/"metadata.json · metadata.csv<br/>all_results_fpmax.csv"/]'); L.append("  minerar_fpmax --> F_meta")
    L.append('  gerar_cql["gerar-cql"]'); L.append('  F_cql[/"neo4j_commands_tree_dataset.cql"/]')
    L.append("  F_meta --> gerar_cql"); L.append("  gerar_cql --> F_cql")
    L.append('  graficos_suporte["graficos-suporte"]'); L.append('  F_pl[/"support_frequence_*.png"/]')
    L.append("  F_meta --> graficos_suporte"); L.append("  graficos_suporte --> F_pl")
    L.append('  consolidar_manifesto["consolidar-manifesto"]')
    L.append('  F_man[/"manifest.json (v2) · manifest.prov.json<br/>log_setup_*.log · output_log.txt"/]')
    for f in ("F_hm", "F_cql", "F_pl"): L.append(f"  {f} --> consolidar_manifesto")
    L.append("  consolidar_manifesto --> F_man")
    L.append("  classDef falha stroke-dasharray: 5 5")
    fal = [mid(n) for n in G if G.nodes[n].get("falha")]
    if fal: L.append("  class " + ",".join(fal) + " falha")
    return "\n".join(L) + "\n"

if __name__ == "__main__" and len(sys.argv) > 2:
    if sys.argv[2] == "--yaml": print(yaml_esboco(sys.argv[1]))
    if sys.argv[2] == "--mermaid": print(mermaid_expandido(sys.argv[1]))
