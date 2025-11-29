#!/usr/bin/env nextflow

params.output_dir = "output"
params.output_file = "hello_world.txt"
params.message = "hello_world"

workflow {
    hello_world()
}

process hello_world {

    // Track the output file
    output:
    path params.output_dir

    // Make it visible in the project folder
    publishDir "${workflow.projectDir}", mode: 'copy'

    script:
    """
    mkdir "${params.output_dir}"
    echo "${params.message}" > "${params.output_dir}/${params.output_file}"
    """
}

// workflow {
//     hello_world()
// }

// process hello_world {

//     // Track the output file
//     output:
//     path "hello.txt"

//     // Make it visible in the project folder
//     publishDir "${workflow.projectDir}", mode: 'copy'

//     script:
//     """
//     echo "hello world" > "hello.txt"
//     """
// }