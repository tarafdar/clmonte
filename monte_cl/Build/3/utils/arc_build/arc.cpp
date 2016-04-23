
#include <stdio.h>
#include <stdlib.h>
#include <cstring>
#include <string>
#include <assert.h>
#include "utils.h"

#define TRUNK_RUN_DIR "bin/host"

// Noted
#define TRUNK_RUN_DIR_ENV "MONTE_RUN=bin/host"
#define TRUNK_BIN_DIR_ENV "MONTE_BIN=bin"
#define TRUNK_SRC_DIR_ENV "MONTE_SRC=monte_cl"
#define TRUNK_TEST_DIR_ENV "MONTE_TEST=test"

// 
#define BUILD_DIR "Build"

// Noted
#define BUILD_RUN_DIR "bin/host"
#define BUILD_BIN_DIR "bin"
#define BUILD_SRC_DIR "monte_cl"
#define BUILD_TEST_DIR "test"

// One for each of the noted above (build or trunk not both)
#define NUM_NEW_ENVS 4

#define SHELL "bash"

void put_env(const char* env)
{
	// Compile errors on putenv if not convert from const char to non const
	char c_path[4096];
	assert(strlen(env) < 4096);
	strcpy(c_path, env);
	c_path[4095] = '\0';
	assert(putenv(c_path) == 0);
}

void add_to_path_env(std::string path_to_add)
{
	char* path = getenv("PATH");
	if(path == NULL) {
		printf("Error: PATH environment variable not set!\n");
		exit(0);
	}

	std::string new_path("PATH=");
	new_path += path;
	new_path += ":";
	new_path += path_to_add;

	put_env(new_path.c_str());
}

int main(int argc, char** argv, char** envp)
{
	printf("Welcome to ARC\n");

	if(argc < 2) {
		printf("Usage:\n\tarc shell\n\tarc commit\n");
		exit(0);
	}

	char* root = getenv("MONTE_ROOT");
	if(root == NULL) {
		printf("Error: MONTE_ROOT environment variable not set!\n");
		exit(0);
	}

	char monte_run[256];
	char monte_bin[256];
	char monte_src[256];
	char monte_test[256];

	if(strcmp(argv[1], "shell") == 0) {
		char* bin = getenv("MONTE_BIN");
		if(bin != NULL) {
			printf("Must exit current shell before entering a new one!\n");
			exit(0);
		}
	
		if(argc != 3) {
			printf("Usage:\n\tarc shell MONTE_KIT\n\tarc shell MONTE_KIT/<build_num>\n");
			exit(0);
		}
		else if(strcmp(argv[2], "MONTE_KIT") == 0) {
			printf("Entering Shell with resource MONTE_KIT (Trunk)\n");
	
			strcpy(monte_run, TRUNK_RUN_DIR_ENV);
			strcpy(monte_bin, TRUNK_BIN_DIR_ENV);
			strcpy(monte_src, TRUNK_SRC_DIR_ENV);
			strcpy(monte_test, TRUNK_TEST_DIR_ENV);

			std::string run_path = std::string(root) + "/" + TRUNK_RUN_DIR;
			add_to_path_env(run_path);
		}
		else if(strncmp(argv[2], "MONTE_KIT/", 10) == 0) {
			char build_num[256];
			sscanf(argv[2] + 10, "%s", build_num);
			if(strlen(build_num) == 0) {
				printf("Error: No build number provided!");
			}
			else {
				printf("Entering Shell with resource MONTE_KIT/%s\n", build_num);
	
				Path build_run_dir_env("MONTE_RUN=");
				build_run_dir_env += root;
				build_run_dir_env.cd(BUILD_DIR);
				build_run_dir_env.cd(build_num);
				build_run_dir_env.cd(BUILD_BIN_DIR);

				Path build_bin_dir_env("MONTE_BIN=");
				build_bin_dir_env += root;
				build_bin_dir_env.cd(BUILD_DIR);
				build_bin_dir_env.cd(build_num);
				build_bin_dir_env.cd(BUILD_BIN_DIR);
	
				Path build_src_dir_env("MONTE_SRC=");
				build_src_dir_env += root;
				build_src_dir_env.cd(BUILD_DIR);
				build_src_dir_env.cd(build_num);
				build_src_dir_env.cd(BUILD_SRC_DIR);

				Path build_test_dir_env("MONTE_TEST=");
				build_test_dir_env += root;
				build_test_dir_env.cd(BUILD_DIR);
				build_test_dir_env.cd(build_num);
				build_test_dir_env.cd(BUILD_TEST_DIR);
				
				Path run_path(root);
				run_path.cd(BUILD_DIR);
				run_path.cd(build_num);
				run_path.cd(BUILD_RUN_DIR);
	
				strcpy(monte_bin, build_bin_dir_env.c_str());
				strcpy(monte_src, build_src_dir_env.c_str());
				strcpy(monte_test, build_test_dir_env.c_str());
				add_to_path_env(std::string(run_path.c_str()));
			}
		}
		else {
			printf("Usage:\n\tarc shell MONTE_KIT\n\tarc shell MONTE_KIT/<build_num>\n");
			exit(0);
		}
	}	
	else if(strcmp(argv[1], "commit") == 0) {
		char* monte_tools = getenv("MONTE_TOOLS");
		if(monte_tools == NULL) {
			printf("Error: MONTE_TOOLS environment variable not set!\n");
			exit(0);
		}

		Path pl_path(monte_tools);
		pl_path.cd("commit.pl");
		std::string cmd = ("perl ");
		cmd += pl_path.to_string();

		system(cmd.c_str());
		exit(0);
	}
	else {
		printf("WTF\n");
		exit(0);
	}

	unsigned int i = 0;
	char** envs = get_expanded_env(envp, NUM_NEW_ENVS, i);

	envs[i++] = monte_run;
	envs[i++] = monte_bin;
	envs[i++] = monte_src;
	envs[i++] = monte_test;
	envs[i] = NULL;

	char* args[2];
	args[0] = SHELL;
	args[1] = NULL;

	execve("/bin/bash", args, envs);

	return 0;
}
