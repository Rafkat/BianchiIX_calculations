.DEFAULT_GOAL := calculate

calculate:
	g++ main.cpp -o ${output_name}
	./${output_name} > ${output_name}.txt
	cp ${output_name}.txt ../ODE_visualisations/{output_name}.txt
