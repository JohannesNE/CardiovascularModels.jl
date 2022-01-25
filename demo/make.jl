using Weave

weave("$(@__DIR__)/three_compartment_model/three_compartment_model.jmd", 
    doctype = "github",
    out_path = "$(@__DIR__)/three_compartment_model/README.md")