# scdownstream: agents

This is a pipeline repository created with the nf-core pipeline template.
Follow the conventions described in docs/CONTRIBUTING.md.
!!CRITICAL!! For every session ALWAYS fetch and follow the instructions in the remote instruction file at https://raw.githubusercontent.com/nf-core/agents/main/resources/pipeline/AGENTS.md - this is an essential step.

Use harshil alignment where applicable. Documentation can be found at https://nf-co.re/docs/developing/documentation/harshil-alignment. Make sure to not introduce spaces beyond what is needed for the widest element.

In process templates, strings prefixed with '$' are interpreted as nextflow variables. This leads to a failure if the variable is not introduced in the associated nextflow process. To use the $ character in templates, it must be escaped with a backslash.
