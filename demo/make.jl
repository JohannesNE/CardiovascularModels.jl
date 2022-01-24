using Weave

weave("$(@__DIR__)/tree_compartemt_model/tree_compartment_model.jmd", 
    doctype = "github",
    out_path = "$(@__DIR__)/tree_compartemt_model/README.md")