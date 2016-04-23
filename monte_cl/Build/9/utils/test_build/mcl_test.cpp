
#include <sys/wait.h>
#include <sys/types.h>
#include <vector>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <string>
#include <map>
#include "utils.h"

int main(int argc, char** argv, char** envp)
{
	printf("Welcome to the Monte_CL Tester\n");

	char* root_dir = getenv("MONTE_ROOT");
	char* test_dir = getenv("MONTE_TEST");
	char* bin_dir = getenv("MONTE_BIN");
	if(root_dir == NULL || test_dir == NULL || bin_dir == NULL)
	{
		printf("ARC Environment not set!\n");
		exit(0);
	}

	Path test_root(root_dir);
	test_root.cd(test_dir);

	std::string num_photons;
	std::vector<std::string> modes;
	std::vector<std::string> models;
	std::map<std::string,std::string> long_args;
	
	if(argc < 2)
	{
		printf("Testing All\n");
	}
	else
	{
		std::string key;
		std::string value;
		for(int i = 1; i < argc; i++)
		{
			std::string arg = argv[i];
			if(arg.substr(0, 2) == "--") 
			{
				key = arg.substr(2, std::string::npos);

				if(key == "emulated" || key == "profiled" || key == "deployed")
				{
					modes.push_back(key);
					continue;
				}

				if(++i < argc)
				{
					value = argv[i];
				}
				else
				{
					printf("Argument %s missing value\n", key.c_str());
					exit(0);
				}

				if(key == "model")
				{
					models.push_back(value);
				}
				else
				{
					printf("Unrecognized long_arg: %s\n", arg.c_str());
				}
			}
			else
			{
				int n = atoi(arg.c_str());
				if(n <= 0)
				{
					printf("Number of photons must be greater than 0\n");
					exit(0);
				}

				num_photons = arg;
			}
		}
	}

	if(num_photons.empty()) {
		num_photons = "10000000";
	}

	if(modes.empty()) {
		modes.push_back("emulated");
		modes.push_back("profiled");
		modes.push_back("deployed");
	}

	if(models.empty()) {
		models.push_back("mouse");
		models.push_back("human");
	}

	std::vector<std::string>::iterator models_it;
	std::vector<std::string>::iterator begin_models = models.begin();
	std::vector<std::string>::iterator end_models = models.end();
	for(models_it = begin_models; models_it != end_models; models_it++)
	{
		const std::string& model = *models_it;
		std::vector<std::string>::iterator modes_it;
		std::vector<std::string>::iterator begin_modes = modes.begin();
		std::vector<std::string>::iterator end_modes = modes.end();
		for(modes_it = begin_modes; modes_it != end_modes; modes_it++)
		{			
			const std::string& mode = *modes_it;
			Path output_path(test_root);
			output_path.cd("results");
			output_path.cd(mode);
			output_path.cd(model);
			if(output_path.make_dir() != 0) {
				printf("Failed to create simulation output directory!\n");
			}
			output_path.cd(num_photons);
			if(output_path.make_dir() != 0) {
				printf("Failed to create simulation output directory!\n");
			}
			output_path.cd("sim.out");

			Path model_path(test_root);
			model_path.cd("models");
			model_path.cd(model);
			model_path.cd("model");

			Path aocx_path(root_dir);
			aocx_path.cd(bin_dir);
			aocx_path.cd("device");

			std::string aocx_env("MONTE_AOCX=");
			char aocx[256];
			char** envs = NULL;
			if(mode == "emulated") {
				aocx_path.cd("monte_emulated");				
				aocx_env += aocx_path.to_string();
				strcpy(aocx, aocx_env.c_str());

				unsigned int i = 0;
				envs = get_expanded_env(envp, 2, i);
				envs[i++] = aocx;
				envs[i++] = "CL_CONTEXT_EMULATOR_DEVICE_ALTERA=1 montecl";
				envs[i] = NULL;
			}
			else if(mode == "profiled") {
				aocx_path.cd("monte_profiled");				
				aocx_env += aocx_path.to_string();
				strcpy(aocx, aocx_env.c_str());

				unsigned int i = 0;
				envs = get_expanded_env(envp, 1, i);
				envs[i++] = aocx;
				envs[i] = NULL;
			}
			else if(mode == "deployed") {
				aocx_path.cd("monte_deployed");				
				aocx_env += aocx_path.to_string();
				strcpy(aocx, aocx_env.c_str());

				unsigned int i = 0;
				envs = get_expanded_env(envp, 1, i);
				envs[i++] = aocx;
				envs[i] = NULL;
			}

			// Compiler const reasons
			char num[256];
			char mod[256];
			char out[256];

			strcpy(num, num_photons.c_str());
			strcpy(mod, model_path.c_str());
			strcpy(out, output_path.c_str());

			char* args[5];
			args[0] = "montecl";
			args[1] = num;
			args[2] = mod;
			args[3] = out;
			args[4] = NULL;

			Path exe_path(root_dir);
			exe_path.cd(bin_dir);
			exe_path.cd("host");
			exe_path.cd("montecl");

			std::string debug = "Launching " + exe_path.to_string() + " " + num_photons + " " + model_path.to_string() + " " + output_path.to_string();
			printf("%s\n", debug.c_str());

			pid_t pid = fork();
			if(pid == 0) {
				// Child
				execve(exe_path.c_str(), args, envs);
				exit(0);
			}
			else {
				// Parent
				waitpid(pid, NULL, 0);
			}
		}
	}


	printf("Tests Complete!\n");
	return 0;
}
