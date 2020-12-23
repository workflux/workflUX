janis inputs test_files/janis/echo_with_graph_changing_input.py -o scratch/test_out/params.yaml

janis translate test_files/janis/echo_with_graph_changing_input.py cwl --inputs test_files/janis/params_string2_enabled.yaml --output-dir scratch/test_out/
janis translate test_files/janis/echo_with_graph_changing_input.py cwl --inputs test_files/janis/params_string2_disabled.yaml --output-dir scratch/test_out/