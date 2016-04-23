
#include <sys/stat.h>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <cerrno>

class Path
{
public:
	Path(std::string str) : str(str) {}
	Path(char* str) : str(str) {}

	void cd(std::string t)
	{
		*this += "/";
		*this += t;
	}

	void cd(const char* t)
	{
		*this += "/";
		*this += t;
	}

	void cd(unsigned int t)
	{
		*this += "/";
		*this += t;
	}

	const std::string& to_string()
	{
		return str;
	}

	const char* c_str()
	{
		return str.c_str();
	}

	int make_dir() const
	{
		int result = mkdir(str.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
		return (result == 0 || errno == EEXIST) ? 0: result;
	}

	template <class T>
	Path& operator+=(T t) {
		str += t;
		return *this;
	}

private:
	std::string str;
};

char** get_expanded_env(char** envp, unsigned int delta_size, unsigned int& orig_size)
{
	char** envp_counter = envp;
	for(orig_size = 0; *envp_counter != NULL; orig_size++) {
		envp_counter++;
	}

	unsigned int new_env_size = orig_size + delta_size;
	char** envs = (char**)malloc(sizeof(char*)*(new_env_size + 1));

	int i;
	for(i = 0; i < orig_size;i++) {
		envs[i] = envp[i];
	}

	envs[new_env_size] = NULL;
	return envs;
}
