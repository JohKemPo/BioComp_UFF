
    MERGE (u:User {uid: '<<USER_UID>>'})
    CREATE (t:Tree {name: 'tree_mafft_iqtree', uid: '<<USER_UID>>', project: 'Projeto_Fixture_G1'})
    MERGE (u)-[:OWNS]->(t);
    

    MATCH (parent {name:'tree_mafft_iqtree', uid: '<<USER_UID>>'})
    WHERE parent:Tree OR parent:Subtree
    CREATE (child:Subtree {name: 'tree_mafft_iqtree_Inner1', uid: '<<USER_UID>>'})
    CREATE (parent)-[:HAS_SUBTREE]->(child);
    

        MATCH (child:Subtree {name: 'tree_mafft_iqtree_Inner1', uid: '<<USER_UID>>'})
        MERGE (s:Support {value: 0.5})
        CREATE (child)-[:HAS_SUPPORT]->(s);
        
MATCH (child:Subtree {name: 'tree_mafft_iqtree_Inner1', uid: '<<USER_UID>>'})
CREATE (m:Metadata {newick_id: "MN000001", terminal_hash: 1111111111111111})
CREATE (child)-[:HAS_METADATA]->(m)
;
MATCH (child:Subtree {name: 'tree_mafft_iqtree_Inner1', uid: '<<USER_UID>>'})
CREATE (m:Metadata {molecule_type: "RNA", organism: "Zika virus", description: "Zika virus isolate X, complete genome", newick_id: "MN000003", terminal_hash: 3333333333333333})
CREATE (child)-[:HAS_METADATA]->(m)
;

    MATCH (parent {name:'tree_mafft_iqtree', uid: '<<USER_UID>>'})
    WHERE parent:Tree OR parent:Subtree
    CREATE (child:Subtree {name: 'tree_mafft_iqtree_Inner2', uid: '<<USER_UID>>'})
    CREATE (parent)-[:HAS_SUBTREE]->(child);
    
MATCH (child:Subtree {name: 'tree_mafft_iqtree_Inner2', uid: '<<USER_UID>>'})
CREATE (m:Metadata {newick_id: "MN000002", terminal_hash: 2222222222222222})
CREATE (child)-[:HAS_METADATA]->(m)
;
MATCH (child:Subtree {name: 'tree_mafft_iqtree_Inner2', uid: '<<USER_UID>>'})
CREATE (m:Metadata {molecule_type: "RNA", organism: "Zika virus", description: "Zika virus isolate X, complete genome", newick_id: "MN000003", terminal_hash: 3333333333333333})
CREATE (child)-[:HAS_METADATA]->(m)
;

    MERGE (u:User {uid: '<<USER_UID>>'})
    CREATE (t:Tree {name: 'tree_mafft_fasttree', uid: '<<USER_UID>>', project: 'Projeto_Fixture_G1'})
    MERGE (u)-[:OWNS]->(t);
    

    MATCH (parent {name:'tree_mafft_fasttree', uid: '<<USER_UID>>'})
    WHERE parent:Tree OR parent:Subtree
    CREATE (child:Subtree {name: 'tree_mafft_fasttree_Inner1', uid: '<<USER_UID>>'})
    CREATE (parent)-[:HAS_SUBTREE]->(child);
    

        MATCH (child:Subtree {name: 'tree_mafft_fasttree_Inner1', uid: '<<USER_UID>>'})
        MERGE (s:Support {value: 0.5})
        CREATE (child)-[:HAS_SUPPORT]->(s);
        

        MATCH (child:Subtree {name: 'tree_mafft_fasttree_Inner1', uid: '<<USER_UID>>'})
        MERGE (s:Support {value: 1.0})
        CREATE (child)-[:HAS_SUPPORT]->(s);
        
MATCH (child:Subtree {name: 'tree_mafft_fasttree_Inner1', uid: '<<USER_UID>>'})
CREATE (m:Metadata {molecule_type: "RNA", topology: "linear", date: "01-JAN-2016", source: "Zika virus", organism: "Zika virus", taxonomy: ["Viruses", "Riboviria", "Flaviviridae", "Orthoflavivirus"], description: "Zika virus isolate d'Ivoire; serogroup: \"Spondweni\", complete genome", newick_id: "MN000001", terminal_hash: 1111111111111111})
CREATE (child)-[:HAS_METADATA]->(m)
CREATE (m)-[:HAS_FEATURE]->(f_0:Feature {type: "source", location: "[0:10807](+)", strand: 1})
CREATE (f_0)-[:HAS_QUALIFIER]->(q_0_0:Qualifier {key: 'organism', value: ["Zika virus"]})
CREATE (f_0)-[:HAS_QUALIFIER]->(q_0_1:Qualifier {key: 'country', value: ["Cote d'Ivoire: Abidjan"]})
CREATE (f_0)-[:HAS_QUALIFIER]->(q_0_2:Qualifier {key: 'note', value: ["passage 1; Vero cells"]})
CREATE (m)-[:HAS_FEATURE]->(f_1:Feature {type: "CDS", location: "[106:10378](+)", strand: 1})
CREATE (f_1)-[:HAS_QUALIFIER]->(q_1_0:Qualifier {key: 'product', value: ["polyprotein"]})
;
MATCH (child:Subtree {name: 'tree_mafft_fasttree_Inner1', uid: '<<USER_UID>>'})
CREATE (m:Metadata {molecule_type: "RNA", topology: "linear", date: "02-FEB-2017", organism: "Zika virus", taxonomy: ["Viruses", "Riboviria"], description: "Zika virus strain ZIKV/Homo_sapiens/BRA//2016, partial genome", newick_id: "MN000002", terminal_hash: 2222222222222222})
CREATE (child)-[:HAS_METADATA]->(m)
CREATE (m)-[:HAS_FEATURE]->(f_0:Feature {type: "source", location: "[0:10000](+)", strand: 1})
CREATE (f_0)-[:HAS_QUALIFIER]->(q_0_0:Qualifier {key: 'country', value: ["Brazil"]})
CREATE (f_0)-[:HAS_QUALIFIER]->(q_0_1:Qualifier {key: 'collection_date', value: ["2016"]})
CREATE (m)-[:HAS_FEATURE]->(f_1:Feature {type: "misc_feature", location: "[5:50]"})
;
MATCH (child:Subtree {name: 'tree_mafft_fasttree_Inner1', uid: '<<USER_UID>>'})
CREATE (m:Metadata {molecule_type: "RNA", organism: "Zika virus", description: "Zika virus isolate X, complete genome", newick_id: "MN000003", terminal_hash: 3333333333333333})
CREATE (child)-[:HAS_METADATA]->(m)
;

    MATCH (parent {name:'tree_mafft_fasttree', uid: '<<USER_UID>>'})
    WHERE parent:Tree OR parent:Subtree
    CREATE (child:Subtree {name: 'tree_mafft_fasttree_Inner2', uid: '<<USER_UID>>'})
    CREATE (parent)-[:HAS_SUBTREE]->(child);
    
MATCH (child:Subtree {name: 'tree_mafft_fasttree_Inner2', uid: '<<USER_UID>>'})
CREATE (m:Metadata {molecule_type: "RNA", topology: "linear", date: "01-JAN-2016", source: "Zika virus", organism: "Zika virus", taxonomy: ["Viruses", "Riboviria", "Flaviviridae", "Orthoflavivirus"], description: "Zika virus isolate d'Ivoire; serogroup: \"Spondweni\", complete genome", newick_id: "MN000001", terminal_hash: 1111111111111111})
CREATE (child)-[:HAS_METADATA]->(m)
CREATE (m)-[:HAS_FEATURE]->(f_0:Feature {type: "source", location: "[0:10807](+)", strand: 1})
CREATE (f_0)-[:HAS_QUALIFIER]->(q_0_0:Qualifier {key: 'organism', value: ["Zika virus"]})
CREATE (f_0)-[:HAS_QUALIFIER]->(q_0_1:Qualifier {key: 'country', value: ["Cote d'Ivoire: Abidjan"]})
CREATE (f_0)-[:HAS_QUALIFIER]->(q_0_2:Qualifier {key: 'note', value: ["passage 1; Vero cells"]})
CREATE (m)-[:HAS_FEATURE]->(f_1:Feature {type: "CDS", location: "[106:10378](+)", strand: 1})
CREATE (f_1)-[:HAS_QUALIFIER]->(q_1_0:Qualifier {key: 'product', value: ["polyprotein"]})
;
MATCH (child:Subtree {name: 'tree_mafft_fasttree_Inner2', uid: '<<USER_UID>>'})
CREATE (m:Metadata {molecule_type: "RNA", topology: "linear", date: "02-FEB-2017", organism: "Zika virus", taxonomy: ["Viruses", "Riboviria"], description: "Zika virus strain ZIKV/Homo_sapiens/BRA//2016, partial genome", newick_id: "MN000002", terminal_hash: 2222222222222222})
CREATE (child)-[:HAS_METADATA]->(m)
CREATE (m)-[:HAS_FEATURE]->(f_0:Feature {type: "source", location: "[0:10000](+)", strand: 1})
CREATE (f_0)-[:HAS_QUALIFIER]->(q_0_0:Qualifier {key: 'country', value: ["Brazil"]})
CREATE (f_0)-[:HAS_QUALIFIER]->(q_0_1:Qualifier {key: 'collection_date', value: ["2016"]})
CREATE (m)-[:HAS_FEATURE]->(f_1:Feature {type: "misc_feature", location: "[5:50]"})
;

    MERGE (u:User {uid: '<<USER_UID>>'})
    CREATE (t:Tree {name: 'tree_mafft_raxml', uid: '<<USER_UID>>', project: 'Projeto_Fixture_G1'})
    MERGE (u)-[:OWNS]->(t);
    

    MATCH (parent {name:'tree_mafft_raxml', uid: '<<USER_UID>>'})
    WHERE parent:Tree OR parent:Subtree
    CREATE (child:Subtree {name: 'tree_mafft_raxml_Inner1', uid: '<<USER_UID>>'})
    CREATE (parent)-[:HAS_SUBTREE]->(child);
    
MATCH (child:Subtree {name: 'tree_mafft_raxml_Inner1', uid: '<<USER_UID>>'})
CREATE (m:Metadata {newick_id: "MN000001", terminal_hash: 1111111111111111})
CREATE (child)-[:HAS_METADATA]->(m)
;
MATCH (child:Subtree {name: 'tree_mafft_raxml_Inner1', uid: '<<USER_UID>>'})
CREATE (m:Metadata {newick_id: "MN000003", terminal_hash: 3333333333333333})
CREATE (child)-[:HAS_METADATA]->(m)
;