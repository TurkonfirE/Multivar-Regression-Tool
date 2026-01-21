#include <iostream>
#include <string>
#include <getopt.h> 
#include "model.h"

using namespace std;

// Samuel Jennings 

int main(int argc, char * argv[]) {
    cout << "Statistical Modeler Program.\n" << std::endl;

    std::ios_base::sync_with_stdio(false);

    bool verbose = false;
    string filename;
    bool file = false;

    option long_options[] = {
        {"help", no_argument, nullptr, 'h'},
        {"verbose", no_argument, nullptr, 'v'},
        {"filename", required_argument, nullptr, 'f'},
        { nullptr, 0, nullptr, '\0' }
    };

    while (1) {
        int option_index = 0;

        int c = getopt_long (argc, argv, "hvf:", long_options, &option_index);

        if (c == -1) break;

        switch (c) {
 
        case 'h':
            cout << "Valid commands include:\n-v (verbose)\n-f *inputfile";
            break;

        case 'v':
            verbose = true;
            break;
        
        case 'f':
            filename = optarg;
            file = true;
        }
    }

    if (file) {
        if (verbose == true) {
            Model model(true, filename); 
            model.make_model();
        }
        else {
            Model model(false, filename); 
            model.make_model();
        }
    }
    else {
        __throw_logic_error("No input file was chosen!");
    }
    
}

