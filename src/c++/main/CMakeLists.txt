# Define all executables with their source files and descriptions
set(HAPPY_EXECUTABLES
    # Core comparison tools
    "xcmp;xcmp.cpp;Compares VCFs by genotypes (truth vs query)"
    "scmp;scmp.cpp;Compares VCFs by alleles (ignoring zygosity)"
    "hapcmp;hapcmp.cpp;Compares VCFs based on sequence alignments"
    
    # Preprocessing tools
    "preprocess;preprocess.cpp;Performs variant decomposition and leftshifting"
    "alleles;alleles.cpp;Pre-processor for scmp that writes one record per allele"
    
    # Enumeration tools
    "hapenum;hapenum.cpp;Takes a VCF and outputs all possible haplotypes"
    "dipenum;dipenum.cpp;Takes a VCF and outputs paired diploid haplotypes"
    
    # Analysis tools
    "quantify;quantify.cpp;Counts variants in VCFs"
    "roc;roc.cpp;Calculates ROC curve from a labeled table with values"
    "blocksplit;blocksplit.cpp ../include/helpers/Roc.hh;Splits VCFs into non-breaking haplotype blocks"
    
    # Utility tools
    "multimerge;multimerge.cpp;Custom version of bcftools merge"
    "vcfhdr2json;vcfhdr2json.cpp;Converts VCF header to JSON"
    "validatevcf;validatevcf.cpp;Checks ref alleles and haplotype enumeration"
    "vcfcheck;vcfcheck.cpp;Tests if VCF is acceptable by hap.py"
    "fastainfo;fastainfo.cpp;Returns contig lengths in JSON format"
    "gvcf2bed;gvcf2bed.cpp;Converts gVCF to BED with regions of homref/variants"
)

# Set common properties for all executables
foreach(exec_info ${HAPPY_EXECUTABLES})
    # Parse the executable information
    string(REPLACE ";" " " exec_info_str "${exec_info}")
    string(REGEX REPLACE "^([^;]+);([^;]+);(.*)$" "\\1;\\2;\\3" exec_parts "${exec_info}")
    
    list(GET exec_parts 0 exec_name)
    list(GET exec_parts 1 exec_sources)
    list(GET exec_parts 2 exec_description)
    
    # Convert space-separated sources to list
    string(REPLACE " " ";" exec_source_list "${exec_sources}")
    
    # Add the executable using our standardized macro
    happy_add_executable(${exec_name} "${exec_source_list}")
    
    # Add target description as a property
    set_target_properties(${exec_name} PROPERTIES 
        DESCRIPTION "${exec_description}"
    )
    
    # Setup standard includes and link with libraries
    happy_setup_includes(${exec_name})
    target_link_libraries(${exec_name} ${HAPLOTYPES_ALL_LIBS})
    
    # Install the executable
    install(TARGETS ${exec_name}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    )
endforeach()
