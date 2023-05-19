#include <string>
#include <unistd.h>
#include <iostream>
#include <sys/wait.h>



int main(int argc, char* argv[]) {

    
    pid_t exeId = fork(); // The child process!

    // If -1, then there was an error in the making of the child!
    if (exeId < 0) {
        std::cerr << "Child process failed!" << std::endl;
    }
    // Child process!
    else if (exeId == 0) {
        //std::cout << "Child here!" << std::endl;

        execl("/bin/ls", "ls", "-l", "-a", NULL);
        return 1;
    }
    // Parent process, stay here!
    else {
        // Wait for child process, then call out that we are in the main program!
        while (wait(NULL) > 0);
        std::cout << "Main program!" << std::endl;
    }
    

}