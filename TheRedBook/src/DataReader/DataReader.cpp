#include "DataReader.h"

ShapeData DataReader::readPLY(std::string path) {

	std::ifstream file(path);
	std::string fileString;
	int vertcount, facecount;
	int vertline = 0, faceline = 0;
	bool header = true;
	ShapeData shapeData;
	bool isPly = false, normalizeColors = false;
	std::vector<std::string> sourceMap, destinationMap;
	std::vector<glm::vec3> vertices, edges, normals;


	// this vector contains all expected values with expected positions in the array
	destinationMap = { {"x"},{"y"},{"z"},{"nx"},{"ny"},{"nz"},{"r"},{"g"},{"b"},{"red"},{"green"},{"blue"},{"a"},{"s"},{"t"} };

	if (!file.good()) {
		std::cout << "file not good" << std::endl;
		std::cout << path << std::endl;
		//ASSERT(false);
		return ShapeData();
	}

	file.seekg(0, std::ios::end);
	fileString.reserve(file.tellg());
	file.seekg(0, std::ios::beg);

	fileString.assign((std::istreambuf_iterator<char>(file)),
		std::istreambuf_iterator<char>());

	std::stringstream ss(fileString);
	std::string to;


	while (std::getline(ss, to, '\n')) {

		if (header) {

			// 0 means "found string"
			//check file for ply string to check if its ply format
			if (!(to.find("ply") == 0) && !isPly) {
				std::cout << "file not in ply format";
				break;
			}
			else if ((to.find("ply") == 0)) {
				isPly = true;
				continue;
			}

			isPly = true;

			// filter lines starting with "format"
			if ((to.find("format") == 0)) {
				continue;
			}

			// filter out lines starting with "comment"
			if (to.find("comment") == 0)
				continue;

			if (to.find("element") == 0) {

				//size for the words vertex / index plus spacebar
				int sizeVertexString = 7;
				int sizeIndexString = 5;

				size_t foundvertex = to.find("vertex");
				size_t foundindex = to.find("face");

				if (foundvertex != std::string::npos) {
					vertcount = stoi(to.substr(foundvertex + sizeVertexString, std::string::npos));
				}

				if (foundindex != std::string::npos) {
					facecount = stoi(to.substr(foundindex + sizeIndexString, std::string::npos));
				}
			}

			if (to.find("property") == 0) {
				int sizeFloatString = 6;

				size_t foundvertexfloat = to.find("float");

				if (foundvertexfloat != std::string::npos) {
					std::string substr = to.substr(foundvertexfloat + sizeFloatString, std::string::npos);

					if (substr == "x" || substr == "y" || (substr == "z")) {
						sourceMap.push_back(substr);
						shapeData.vertMapArr[0][0]++;
					}
					if (substr == "nx" || substr == "ny" || (substr == "nz")) {
						sourceMap.push_back(substr);
						shapeData.vertMapArr[1][0]++;
					}
					if (substr == "r" || substr == "g" || substr == "b" || substr == "red" || substr == "green" || substr == "blue" || (substr == "a")) {
						sourceMap.push_back(substr);
						shapeData.vertMapArr[2][0]++;
					}
					if (substr == "s" || substr == "t" || (substr == "u")) {
						sourceMap.push_back(substr);
						shapeData.vertMapArr[3][0]++;
					}
				}

				size_t foundvertexuchar = to.find("uchar");

				if (foundvertexuchar != std::string::npos) {
					std::string substr = to.substr(foundvertexuchar + sizeFloatString, std::string::npos);

					normalizeColors = true;

					if (substr == "x" || substr == "y" || (substr == "z")) {
						sourceMap.push_back(substr);
						shapeData.vertMapArr[0][0]++;
					}
					if (substr == "nx" || substr == "ny" || (substr == "nz")) {
						sourceMap.push_back(substr);
						shapeData.vertMapArr[1][0]++;
					}
					if (substr == "r" || substr == "g" || substr == "b" || substr == "red" || substr == "green" || substr == "blue" || (substr == "a")) {
						sourceMap.push_back(substr);
						shapeData.vertMapArr[2][0]++;
					}
					if (substr == "s" || substr == "t" || (substr == "u")) {
						sourceMap.push_back(substr);
						shapeData.vertMapArr[3][0]++;
					}
				}
			}
		}
		if (to.find("end_header") == 0) {
			header = false;

			shapeData.mapSize = (shapeData.vertMapArr[0][0] + shapeData.vertMapArr[1][0] + shapeData.vertMapArr[2][0] + shapeData.vertMapArr[3][0]);

			continue;
		}
		if (!header) {

			if (vertline < vertcount) {

				std::stringstream ssin(to);

				int offset = 0;
				std::vector<GLfloat> coords, normals, textures, colors;

				for (int i = 0; i < shapeData.mapSize; i++) {

					std::string vert;

					bool found = false;

					ssin >> vert;

					//map values from source file order to the required order
					if (sourceMap[i] == "x" || sourceMap[i] == "y" || sourceMap[i] == "z") {
						coords.push_back(atof(vert.c_str()));
					}
					if (sourceMap[i] == "nx" || sourceMap[i] == "ny" || (sourceMap[i] == "nz")) {
						normals.push_back(atof(vert.c_str()));
					}
					if (sourceMap[i] == "r" || sourceMap[i] == "g" || sourceMap[i] == "b" || sourceMap[i] == "red" || sourceMap[i] == "green" || sourceMap[i] == "blue" || (sourceMap[i] == "a")) {
						colors.push_back(atof(vert.c_str()) / 255);
					}
					if (sourceMap[i] == "s" || sourceMap[i] == "t" || (sourceMap[i] == "u")) {
						textures.push_back(atof(vert.c_str()));
					}

				}

				/*	if (shapeData.vertMapArr[2][0] < 2) {
						for (int i = 0; i < 3; i++)
							colors.push_back(0.0f);
					}*/

				for (GLfloat _coords : coords) {
					shapeData.vertices.push_back(_coords);					
					//std::cout << _coords << " ";
				}

				shapeData.geometryPoints.push_back(glm::vec3(coords[0], coords[1], coords[2]));

				for (GLfloat _normals : normals) {
					shapeData.vertices.push_back(_normals);
					//std::cout << _normals << " ";
				}

				for (GLfloat _colors : colors) {
					shapeData.vertices.push_back(_colors);
					//std::cout << _colors << " ";
				}

				for (GLfloat _textures : textures) {
					shapeData.vertices.push_back(_textures);
					//std::cout << _textures << " ";
				}

				//std::cout << std::endl;

				vertline++;
			}
			else {
				if (faceline < facecount) {

					std::stringstream ssin(to);

					for (int i = 0; i < 4; i++) {

						std::string face;

						ssin >> face;

						//skip the linesize indicator
						if (i == 0)
							continue;

						shapeData.indices.push_back(atoi(face.c_str()));

					}
					faceline++;
				}

			}

		}
	}

	unsigned int offset = 2; // offset for actual vertice and next vertices
	const int edgeCount = shapeData.geometryPoints.size() / offset - 1;
	for (int i = 0; i < (shapeData.geometryPoints.size() / offset - 1); i++) {
		shapeData.edges.push_back(Mathematics::CalcEdge(shapeData.geometryPoints[i * offset], shapeData.geometryPoints[i * offset + 1]));
	}

	offset = shapeData.mapSize;
	for (int i = 0; i < vertcount; i++)
		shapeData.normals.push_back(glm::vec3(shapeData.vertices[i * offset + 3], shapeData.vertices[i * offset + 1 + 3], shapeData.vertices[i * offset + 2 + 3]));

	shapeData.edgeSize = shapeData.edges.size() * sizeof(GLfloat);
	shapeData.edgeCount = shapeData.edges.size();
	shapeData.normalSize = shapeData.normals.size() * sizeof(GLfloat);
	shapeData.normalCount = shapeData.normals.size();

	shapeData.vertexSize = sizeof(GLfloat) * shapeData.mapSize * vertcount;
	shapeData.vertexCount = vertcount;

	shapeData.indexSize = sizeof(GLuint) * facecount * 3;
	shapeData.indexCount = facecount;

	shapeData.boundingBox = Mathematics::GenerateBoundingBox(shapeData.vertices);

	return shapeData;
}