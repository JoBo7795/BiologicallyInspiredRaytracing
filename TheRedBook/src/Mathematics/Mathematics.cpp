#include "Mathematics.h"

BoundBox Mathematics::GenerateBoundingBox(std::vector<float> coordList) {

	BoundBox boundingBox;
	int offset = 5;

	boundingBox.xmin = coordList[0];
	boundingBox.xmax = coordList[0];
	boundingBox.ymin = coordList[1];
	boundingBox.ymax = coordList[1];
	boundingBox.zmin = coordList[2];
	boundingBox.zmax = coordList[2];

	for (int i = 0; i < coordList.size() / offset; i++)
	{
		if (boundingBox.xmax < coordList[i * offset])
			boundingBox.xmax = coordList[i * offset];
		if (boundingBox.xmin > coordList[i * offset])
			boundingBox.xmin = coordList[i * offset];


		if (boundingBox.ymax < coordList[i * offset + 1])
			boundingBox.ymax = coordList[i * offset + 1];
		if (boundingBox.ymin > coordList[i * offset + 1])
			boundingBox.ymin = coordList[i * offset + 1];


		if (boundingBox.zmax < coordList[i * offset + 2])
			boundingBox.zmax = coordList[i * offset + 2];
		if (boundingBox.zmin > coordList[i * offset + 2])
			boundingBox.zmin = coordList[i * offset + 2];
	}

	boundingBox.xList.push_back(boundingBox.xmax);
	boundingBox.xList.push_back(boundingBox.xmin);
	boundingBox.yList.push_back(boundingBox.ymax);
	boundingBox.yList.push_back(boundingBox.ymin);
	boundingBox.zList.push_back(boundingBox.zmax);
	boundingBox.zList.push_back(boundingBox.zmin);

	return boundingBox;
}
#include <iostream>
BoundBox Mathematics::GenerateBoundingBox(glm::vec3 centerPoint, float dimension) {

	BoundBox boundingBox;

	//std::cout << "centerPointX: " << centerPoint.x << " centerPointY: " << centerPoint.y << " centerPointZ: " << centerPoint.z << std::endl;

	boundingBox.xList.reserve(2);
	boundingBox.yList.reserve(2);
	boundingBox.zList.reserve(2);

	boundingBox.xmin = centerPoint.x ;
	boundingBox.xmax = centerPoint.x + dimension;
	boundingBox.ymin = centerPoint.y;
	boundingBox.ymax = centerPoint.y + dimension;
	boundingBox.zmin = centerPoint.z;
	boundingBox.zmax = centerPoint.z + dimension;

	boundingBox.xList.push_back(centerPoint.x - dimension);
	boundingBox.xList.push_back(centerPoint.x + dimension);
	boundingBox.yList.push_back(centerPoint.y - dimension);
	boundingBox.yList.push_back(centerPoint.y - dimension);
	boundingBox.zList.push_back(centerPoint.z - dimension);
	boundingBox.zList.push_back(centerPoint.z + dimension);

	return boundingBox;
}

BoundBox Mathematics::GenerateBoundingBox(std::vector<glm::vec4> coordList) {

	BoundBox boundingBox;
	int offset = 5;

	boundingBox.xmin = coordList[0].x;
	boundingBox.xmax = coordList[0].x;
	boundingBox.ymin = coordList[0].y;
	boundingBox.ymax = coordList[0].y;
	boundingBox.zmin = coordList[0].z;
	boundingBox.zmax = coordList[0].z;

	int size = coordList.size();

	for (int i = 0; i < size; i++)
	{
		if (boundingBox.xmax < coordList[i].x)
			boundingBox.xmax = coordList[i].x;
		if (boundingBox.xmin > coordList[i].x)
			boundingBox.xmin = coordList[i].x;


		if (boundingBox.ymax < coordList[i].y)
			boundingBox.ymax = coordList[i].y;
		if (boundingBox.ymin > coordList[i].y)
			boundingBox.ymin = coordList[i].y;


		if (boundingBox.zmax < coordList[i].z)
			boundingBox.zmax = coordList[i].z;
		if (boundingBox.zmin > coordList[i].z)
			boundingBox.zmin = coordList[i].z;

	}


	boundingBox.xList.reserve(2);
	boundingBox.xList.push_back(boundingBox.xmax);
	boundingBox.xList.push_back(boundingBox.xmin);

	boundingBox.yList.reserve(2);
	boundingBox.yList.push_back(boundingBox.ymax);
	boundingBox.yList.push_back(boundingBox.ymin);

	boundingBox.zList.reserve(2);
	boundingBox.zList.push_back(boundingBox.zmax);
	boundingBox.zList.push_back(boundingBox.zmin);

	return boundingBox;
}

BoundBox Mathematics::GenerateBoundingBox(std::vector<glm::vec3> coordList) {

	BoundBox boundingBox;
	int offset = 5;

	boundingBox.xmin = coordList[0].x;
	boundingBox.xmax = coordList[0].x;
	boundingBox.ymin = coordList[0].y;
	boundingBox.ymax = coordList[0].y;
	boundingBox.zmin = coordList[0].z;
	boundingBox.zmax = coordList[0].z;

	for (int i = 0; i < coordList.size(); i++)
	{
		if (boundingBox.xmax < coordList[i].x)
			boundingBox.xmax = coordList[i].x;
		if (boundingBox.xmin > coordList[i].x)
			boundingBox.xmin = coordList[i].x;


		if (boundingBox.ymax < coordList[i].y)
			boundingBox.ymax = coordList[i].y;
		if (boundingBox.ymin > coordList[i].y)
			boundingBox.ymin = coordList[i].y;


		if (boundingBox.zmax < coordList[i].z)
			boundingBox.zmax = coordList[i].z;
		if (boundingBox.zmin > coordList[i].z)
			boundingBox.zmin = coordList[i].z;

	}

	boundingBox.xList.push_back(boundingBox.xmax);
	boundingBox.xList.push_back(boundingBox.xmin);
	boundingBox.yList.push_back(boundingBox.ymax);
	boundingBox.yList.push_back(boundingBox.ymin);
	boundingBox.zList.push_back(boundingBox.zmax);
	boundingBox.zList.push_back(boundingBox.zmin);

	return boundingBox;
}

glm::vec3 Mathematics::CalcEdge(glm::vec3 vectorA, glm::vec3 vectorB) {


	float edgeVecX = vectorB.x - vectorA.x;
	float edgeVecY = vectorB.y - vectorA.y;
	float edgeVecZ = vectorB.z - vectorA.z;

	return glm::vec3(edgeVecX, edgeVecY, edgeVecZ);
}

void Mathematics::GenerateNormalsFromBB(BoundBox* bb) {

	std::vector<glm::vec3> edges;
	std::vector<glm::vec3>* normals;	

	edges.push_back(glm::normalize(Mathematics::CalcEdge(glm::vec3(bb->xmin, bb->ymin, bb->zmax),
		glm::vec3(glm::vec3(bb->xmax, bb->ymin, bb->zmax)))));

	edges.push_back(glm::normalize(Mathematics::CalcEdge(glm::vec3(bb->xmax, bb->ymin, bb->zmax),
		glm::vec3(glm::vec3(bb->xmax, bb->ymax, bb->zmax)))));
	//left + right

	edges.push_back(glm::normalize(Mathematics::CalcEdge(glm::vec3(bb->xmin, bb->ymin, bb->zmin),
		glm::vec3(glm::vec3(bb->xmin, bb->ymin, bb->zmax)))));

	edges.push_back(glm::normalize(Mathematics::CalcEdge(glm::vec3(bb->xmin, bb->ymin, bb->zmax),
		glm::vec3(glm::vec3(bb->xmin, bb->ymax, bb->zmax)))));
	// up + down

	edges.push_back(glm::normalize(Mathematics::CalcEdge(glm::vec3(bb->xmin, bb->ymin, bb->zmin),
		glm::vec3(glm::vec3(bb->xmax, bb->ymin, bb->zmin)))));
										
	edges.push_back(glm::normalize(Mathematics::CalcEdge(glm::vec3(bb->xmax, bb->ymin, bb->zmin),
		glm::vec3(glm::vec3(bb->xmax, bb->ymin, bb->zmax)))));

	bb->normals[0] = glm::normalize(glm::cross(edges[0], edges[1])); // vorne
	bb->normals[1] = glm::normalize(-1.0f * glm::cross(edges[0], edges[1])); // hinten
	bb->normals[2] = glm::normalize(glm::cross(edges[2], edges[3])); // links
	bb->normals[3] = glm::normalize(-1.0f * glm::cross(edges[2], edges[3])); // rechts
	bb->normals[4] = glm::normalize(glm::cross(edges[4], edges[5])); // unten
	bb->normals[5] = glm::normalize(-1.0f * glm::cross(edges[4], edges[5])); //oben					

	//return normals; // returned leeren vector !!
}

std::vector<Mathematics::KDNode> Mathematics::CreateKDTree(std::vector<std::array<glm::vec3,3>>& primitiveList,int depth, int& root) {

	int DIM = 3;		
	std::vector<glm::vec3> dimMiddleList;
	for (int i = 0; i < primitiveList.size(); i++)
	{
		/*glm::vec3 middleVec;

		middleVec.x = (primitiveList[i][0].x + primitiveList[i][1].x + primitiveList[i][2].x)/3;
		middleVec.y = (primitiveList[i][0].y + primitiveList[i][1].y + primitiveList[i][2].y)/3;
		middleVec.z = (primitiveList[i][0].z + primitiveList[i][1].z + primitiveList[i][2].z)/3;*/

		dimMiddleList.push_back(primitiveList[i][0]);
		dimMiddleList.push_back(primitiveList[i][1]);
		dimMiddleList.push_back(primitiveList[i][2]);

	}

	std::vector<Mathematics::KDNode> nodes;

	// convert points in pointList into node format
	for (int i = 0; i < dimMiddleList.size(); i++) {
		Mathematics::KDNode node;
		node.data = dimMiddleList[i];

		nodes.push_back(node);
	}

	// create tree, return rootIndex
	root = createTree(nodes,0,nodes.size(),0,DIM,depth);

	return nodes;	

}
std::vector<int> primArr;
void Mathematics::SortPrimitivesToNodes(std::vector<Mathematics::KDNode>& nodes, std::vector<std::array<glm::vec3, 3>>& primitiveList,int root, int depth) {
	int dimensions = 3;
	int primCounter = -1;
	for (auto& primitive : primitiveList) {

		int currentNode = root;
		primCounter++;
		for (int i = 0;i< primitive.size(); i++) {
			for (int j = 0; j < depth; j++) {
				if ((nodes[currentNode].data[j % dimensions] < primitive[i][j % dimensions]) && (nodes[currentNode].indexLeft > 0)) {
					currentNode = nodes[currentNode].indexLeft;
				}
				else if ((nodes[currentNode].data[j % dimensions] > primitive[i][j % dimensions]) && (nodes[currentNode].indexRight > 0)) {
					currentNode = nodes[currentNode].indexRight;
				}
				else {
					nodes[currentNode].PrimitivePointerArr.push_back(primCounter * primitive.size() + i);
					nodes[currentNode].FindMinMax(primitive[i]);
					primArr.push_back(primCounter* primitive.size() + i);
					currentNode = root;
					
					break;
				}
			}
		}

	}

}

int Mathematics::createTree(std::vector<Mathematics::KDNode>& nodes,int startIndex, int endIndex, int index, int dimensions, int depth)
{	
	if (startIndex >= endIndex || depth <=0)
	{
		return -1;
	}

	int meanIndex = startIndex + (endIndex - startIndex) / 2; 
	std::vector<Mathematics::KDNode>::iterator i = nodes.begin();
	std::nth_element(i + startIndex, i + meanIndex, i + endIndex, [&](auto a, auto b) {return a.data[index % dimensions] < b.data[index % dimensions]; });
	index = (index + 1) % dimensions;

	nodes[meanIndex].indexLeft = createTree(nodes,startIndex, meanIndex, index, dimensions, depth-1);
	if(nodes[meanIndex].indexLeft > -1)
		nodes[nodes[meanIndex].indexLeft].indexParent = meanIndex;
	nodes[meanIndex].indexRight = createTree(nodes,meanIndex + 1, endIndex, index, dimensions,depth-1);
	if (nodes[meanIndex].indexRight > -1)
		nodes[nodes[meanIndex].indexRight].indexParent = meanIndex;

	return meanIndex;
}



std::vector<float> Mathematics::NodeToFloatList(std::vector<Mathematics::KDNode>& nodeList) {

	std::vector<float> fNodeList;

	for (int i = 0; i < nodeList.size(); i++) {
		fNodeList.push_back(static_cast<float>(nodeList[i].data.x));
		fNodeList.push_back(static_cast<float>(nodeList[i].data.y));
		fNodeList.push_back(static_cast<float>(nodeList[i].data.z));
		fNodeList.push_back(0.0f);

		fNodeList.push_back(static_cast<float>(nodeList[i].indexLeft));
		fNodeList.push_back(static_cast<float>(nodeList[i].indexRight));
		fNodeList.push_back(static_cast<float>(nodeList[i].indexParent));
		fNodeList.push_back(0.0f);

		fNodeList.push_back(static_cast<float>(nodeList[i].minX));
		fNodeList.push_back(static_cast<float>(nodeList[i].minY));
		fNodeList.push_back(static_cast<float>(nodeList[i].minZ));
		fNodeList.push_back(0.0f);

		fNodeList.push_back(nodeList[i].maxX);
		fNodeList.push_back(nodeList[i].maxY);
		fNodeList.push_back(nodeList[i].maxZ);
		fNodeList.push_back(0.0f);

		fNodeList.push_back(static_cast<float>(nodeList[i].PrimRangeL));
		fNodeList.push_back(static_cast<float>(nodeList[i].PrimRangeH));
		fNodeList.push_back(0.0f);
		fNodeList.push_back(0.0f);
	}

	return fNodeList;
}

std::vector<int> Mathematics::CreatePrimitiveArrayList(std::vector<Mathematics::KDNode>& nodeList) {

	std::vector<int> primList;
	int currentIndex = 0, startIndex = 0;
	for (Mathematics::KDNode& node : nodeList) {
		int size = node.PrimitivePointerArr.size();
		if (size > 0) {
			
			for (int i = 0; i < size; i++) {
				primList.push_back(node.PrimitivePointerArr[i]);
				currentIndex++;
			}

			node.PrimRangeL = startIndex;
			node.PrimRangeH = currentIndex;
			startIndex = currentIndex;
		}
		
	}

	return primList;

}
