#include <iostream> 
#include <fstream>
#include <sstream>
#include <string>
#include <experimental/filesystem>
#include <eigen3/Eigen/Dense>

using namespace Eigen;
namespace fs = std::experimental::filesystem;

VectorXd load_pgm(const std::string &filename) {
	// returns a picture as a flattened vector

	int row = 0, col = 0, rows = 0, cols = 0;

	std::ifstream infile(filename);
	std::stringstream ss;
	std::string inputLine = "";

	// First line : version
	std::getline(infile,inputLine);

	// Second line : comment
	std::getline(infile,inputLine);

	// Continue with a stringstream
	ss << infile.rdbuf();
	// Third line : size
	ss >> cols >> rows;

	VectorXd picture(rows*cols);

	// Following lines : data
	for(row = 0; row < rows; ++row) {
		for (col = 0; col < cols; ++col) {
			int val;
			ss >> val;
			picture(col*rows + row) = val;
		}
	}

	infile.close();
	return picture;
}

void recognize(std::string testPicPath) {
	int h = 231;
	int w = 195;
	int M = 15;

	MatrixXd faces(h*w, M);
	VectorXd meanFace(h*w);
	
	// loads pictures as flattened vectors into faces
	for (int i=0; i<M; i++) {
		std::string filename = "./basePictures/subject"+ 
			std::to_string(i+1) + ".pgm";
		VectorXd flatPic = load_pgm(filename);
		faces.col(i) = flatPic;
		// TODO: Point (b)
		meanFace += flatPic;
	}
	// calculate mean face and sub it from all faces column wise
	meanFace /= M;
	faces.colwise() -= meanFace;

	
	// TODO: Point (e)
	JacobiSVD<MatrixXd> svd(faces, ComputeThinU | ComputeThinV);
	MatrixXd U = svd.matrixU(); // This is U but also AAt

	// try to recognize a test face
	std::string testPicName = testPicPath;
	VectorXd newFace = load_pgm(testPicName);

	// TODO: Point (f)
	VectorXd projNewFace = U.transpose() * (newFace-meanFace);
	
	// TODO: Point (g)
	MatrixXd projFaces = U.transpose() * faces;

	int indexMinNorm;
	(projFaces.colwise()-projNewFace).colwise().norm().minCoeff(&indexMinNorm);
	std::cout << testPicName << " is identified as subject "
		<< indexMinNorm+1 << std::endl;

	std::string comm = "gnome-terminal -x sh -c 'eog ./basePictures/subject"+std::to_string(indexMinNorm+1)+".pgm	'";

	system(comm.c_str());
}

int main() {

	bool exit = false;
	std::string path = "./testPictures";
	std::vector<fs::directory_entry> directory;
	for (const auto & entry : fs::directory_iterator(path)) {
		directory.push_back(entry);
	}
	while(!exit){
		std::cout << "--- Choose Test Face ---" << std::endl;
		int count = 0;
		for (const auto & entry : directory) {
			std::cout << count << ". " << entry.path() << std::endl;
			count++;
		}
		int chosen = 0;
		std::cin >> chosen;
		
		if (chosen > directory.size() || chosen < 0){
			exit = true;
			continue;
		}

		std::string chosenPath = directory[chosen].path();

		std::string comm = "gnome-terminal -x sh -c 'eog "+chosenPath+"'";

		system(comm.c_str());
		
		recognize(chosenPath);//"./testPictures/Narutowicz.pgm"
		std::cout << "\n\n" << std::endl;
	}
}
