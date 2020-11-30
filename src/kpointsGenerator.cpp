//# include <kpointsGenerator.h>
# include <Eigen/Dense>
# include <string> 
# include <fstream>
# include <iostream>

using namespace Eigen;
using namespace std;

enum class GeneratorMode:char{
    ReadFromFile,
    GenerateFromFile
};

class kpointsGenerator
{
private:
    Array<double, Dynamic, Dynamic> kpointsArray();

public:
    kpointsGenerator(GeneratorMode ReadFromFile, const string* FolderName) {
        int loopIndex = 1;
        int startReadLine = 4;
        string FileName ("\\IBZKPT");
        string LineContent;
        string Location = *FolderName + FileName;
        fstream instream(Location, ios_base::in);
        if (instream) {
            while (instream)
            {
                if (loopIndex >= startReadLine) break;
                // getline(instream, NULL, '\n');
                loopIndex++;
            };
            
            while (instream)
            {
                getline(instream, LineContent, '\n');
                cout << LineContent;
            }
            
        }
    };
    kpointsGenerator() {
        cout << "constructor Called" << endl;
    };
    ~kpointsGenerator() {
        cout << "Destructor called" << endl;
    };
};


// kpointsGenerator::kpointsGenerator(/* args */)
// {
// }

// kpointsGenerator::~kpointsGenerator()
// {
// }

int main(void) {
    string Location ("D:\\OneDrive - tongji.edu.cn\\vscode_workspace\\DESKTOP-DQVLUVG\\VScode_WorkSpace\\model\\src");
    kpointsGenerator *kpoints;
    GeneratorMode mode = GeneratorMode::ReadFromFile;
    kpoints = new kpointsGenerator(mode, &Location);
    delete kpoints;
    return 0;
}
