#include "DataWriter.h"

void DataWriter::WriteToProgramPathFile() {

	time_t rawtime;
	struct tm timeinfo;
	char dateTime[80];

	time(&rawtime);
	localtime_s(&timeinfo, &rawtime);

	char date[20];
	strftime(date, sizeof(date), "%Y-%m-%d", &timeinfo);

	strftime(dateTime, sizeof(dateTime), "%Y_%m_%d_%H_%M_%S", &timeinfo);
	std::string str(dateTime);
	ProgramParams::path += (date + std::string("\\"));
	std::string num = std::string(dateTime);
	std::string path2 = (char*)".png";
	std::string filepath = ProgramParams::path + ProgramParams::paramPaths + num + ProgramParams::params + path2;
	const char* finalPath = filepath.c_str();
	std::cout << "The image is saved under path: " << finalPath << std::endl;

	if(ProgramParams::CPURender){
		//stbi_write_jpg("C:\\Users\\Johannes\\Desktop\\Renderings\\raytrace24.jpg", ProgramParams::windowWidth, ProgramParams::windowHeight, CHANNEL_NUM, CPUMode::pixels, 100);
		stbi_write_jpg(finalPath, ProgramParams::windowWidth, ProgramParams::windowHeight, CHANNEL_NUM, CPUMode::pixels, 100);
		//ProgramParams::window.saveImageCPU((char*)finalPath, CPUMode::pixels);
	}
	else {
		ProgramParams::window.saveImage((char*)finalPath, ProgramParams::window.GetWindowRef());
	}
}