#include "hlp.h"

#include <fstream>
#include <deal.II/numerics/vector_tools.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <string>
#include <sstream>

using namespace std;


#ifdef _WIN32
    #include <direct.h>
    #define MKDIR(path) _mkdir(path)
#else
    #define MKDIR(path) mkdir(path, 0777)
#endif

bool aspect::DizHlp::ensureFolderExists ( const std::string & folderPath ) {
    struct stat info;

    std::stringstream pathStream(folderPath);
    std::string partialPath;
    std::string segment;

    // Handle absolute paths on Unix/Windows
#ifdef _WIN32
    if (folderPath.size() > 1 && folderPath[1] == ':') {
        partialPath = folderPath.substr(0, 2); // Keep drive letter
        pathStream.seekg(2); // Skip drive letter
    }
#else
    if (folderPath[0] == '/') {
        partialPath = "/";
    }
#endif

    // Parse and create each directory segment
    while (std::getline(pathStream, segment, '/')) {
        if (segment.empty()) continue;

        partialPath += (partialPath.empty() || partialPath.back() == '/') ? segment : ("/" + segment);

        if (stat(partialPath.c_str(), &info) != 0) { 
            // Directory does not exist; create it
            if (MKDIR(partialPath.c_str()) != 0) {
                perror(("Failed to create folder: " + partialPath).c_str());
                return false;
            }
        } else if (!(info.st_mode & S_IFDIR)) {
            // Path exists, but it's not a directory
            std::cerr << partialPath << " exists but is not a directory.\n";
            return false;
        }
    }

    return true;
}
