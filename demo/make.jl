using Weave

weave("$(@__DIR__)/three_compartment_model/three_compartment_model.jmd", 
    doctype = "github",
    out_path = "$(@__DIR__)/three_compartment_model/README.md")

weave("$(@__DIR__)/simple_lung/simple_lung.jmd", 
    doctype = "github",
    out_path = "$(@__DIR__)/simple_lung/README.md")

weave("$(@__DIR__)/parameter_estimation/parameter_estimation.jmd", 
    doctype = "github",
    out_path = "$(@__DIR__)/parameter_estimation/README.md")
    