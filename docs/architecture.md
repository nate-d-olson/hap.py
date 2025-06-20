# Package Architecture

```dot
// Graphviz diagram showing high level components
// Run `dot -Tpng architecture.md -o architecture.png` to render

digraph architecture {
    rankdir=LR;
    node [shape=box, style=filled, fillcolor=lightgray];
    subgraph cluster_cli {
        label="Command line tools";
        hap[label="hap.py"];
        pre[label="pre.py"];
        qfy[label="qfy.py"];
    }
    subgraph cluster_haplo {
        label="hap_py.haplo";
        preprocess[label="python_preprocess"];
        quantify[label="python_quantify"];
        cmp[label="python_hapcmp"];
        vcfcheck[label="python_vcfcheck"];
        vcfeval[label="vcfeval"];
    }
    tools[label="hap_py.tools"];
    utils[label="hap_py.utils"];

    hap -> preprocess;
    hap -> cmp;
    hap -> vcfeval;
    pre -> preprocess;
    qfy -> quantify;

    preprocess -> tools;
    quantify -> tools;
    cmp -> tools;
    vcfeval -> tools;
    vcfeval -> utils;
}
```

This diagram illustrates how the top-level command line interfaces rely on the core modules under `hap_py.haplo`, which in turn use helper functions in `hap_py.tools` and `hap_py.utils`.
