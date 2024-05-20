#pragma once

#include <vector>
#include <algorithm>
#include <array>
#include "../Objects/BBStruct.h"
#include "../../Dependencies/glm/glm/glm.hpp"
#include "../../Dependencies/glm/glm/gtc/matrix_transform.hpp"
#include "../../Dependencies/glm/glm/gtc/type_ptr.hpp"



namespace Mathematics
{
	struct KDNode {
		glm::vec3 data;
		int indexLeft, indexRight, indexParent;
		int PrimRangeL = 0, PrimRangeH = 0;
		float minX, minY, minZ, maxX, maxY, maxZ;

		bool initMinX = false, initMinY = false, initMinZ = false, initMaxX = false, initMaxY = false, initMaxZ = false;
		std::vector<int> PrimitivePointerArr;

		void FindMinMax(glm::vec3 point) {

			if ((point.x < minX) || !initMinX) {
				minX = point.x;
				initMinX = true;
			}
			if ((point.y < minY) || !initMinY) {
				minY = point.y;
				initMinY = true;
			}
			if ((point.z < minZ) || !initMinZ) {
				minZ = point.z;
				initMinZ = true;
			}

			if ((point.x > maxX) || !initMaxX) {
				maxX = point.x;
				initMaxX = true;
			}
			if ((point.y > maxY) || !initMaxY) {
				maxY = point.y;
				initMaxY = true;
			}
			if ((point.z > maxZ) || !initMaxZ) {
				maxZ = point.z;
				initMaxZ = true;
			}
		}
	};

	struct LeafData {
		int NumLeafs;
		//int PrimitivePointerArr[];
		int PrimRangeL, PrimRangeH;
	};

	KDNode Insert(glm::vec3 x, KDNode& t, int cd, int DIM);

	BoundBox GenerateBoundingBox(std::vector<float> coordList);
	BoundBox GenerateBoundingBox(glm::vec3 centerPoint, float dimension);
	BoundBox GenerateBoundingBox(std::vector<glm::vec4> coordList);
	BoundBox GenerateBoundingBox(std::vector<glm::vec3> coordList);
	std::vector<KDNode> CreateKDTree(std::vector<std::array<glm::vec3, 3>>& primitiveList, int depth, int& root);
	int createTree(std::vector<Mathematics::KDNode>& nodes, int startIndex, int endIndex, int index, int dimensions, int depth);
	void SortPrimitivesToNodes(std::vector<Mathematics::KDNode>& nodes, std::vector<std::array<glm::vec3, 3>>& primitiveList, int root,int depth);
	std::vector<float> NodeToFloatList(std::vector<Mathematics::KDNode>& nodeList);
	std::vector<int> CreatePrimitiveArrayList(std::vector<Mathematics::KDNode>& nodeList);

	glm::vec3 CalcEdge(glm::vec3 vectorA, glm::vec3 vectorB);
	void GenerateNormalsFromBB(BoundBox* bb);
};

