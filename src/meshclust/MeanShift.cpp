/*
 MeShClust v3.0 clusters sequences using the mean shift algorithm and alignment-free identity scores.

 Copyright (C) 2020-2022 Hani Z. Girgis, PhD

 Academic use: Affero General Public License version 1.

 Any restrictions to use for-profit or non-academics: Alternative commercial license is needed.

 This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

 Please contact Dr. Hani Z. Girgis (hzgirgis@buffalo.edu) if you need more information.
 */

/*
 * MeanShift.cpp
 *
 *  Created on: Dec 23, 2020
 *      Author: Hani Zakaria Girgis, PhD
 */

template<class V>
MeanShift<V>::MeanShift(
		std::tuple<V**, uint64_t**, std::string**, int*, int> tup,
		IdentityCalculator<V> &id, int t, double th) :
		identity(id), threadNum(t), threshold(th) {

	kHistList = std::get < 0 > (tup);
	monoHistList = std::get < 1 > (tup);
	infoList = std::get < 2 > (tup);
	lenList = std::get < 3 > (tup);
	size = std::get < 4 > (tup);
	isDataCleared = false;
	start();
}

template<class V>
MeanShift<V>::MeanShift(Block *block, IdentityCalculator<V> &id, int t,
		double th) :
		identity(id), threadNum(t), threshold(th) {

	initData(block);
	start();
}

template<class V>
void MeanShift<V>::clearData() {
	identity.freeBlock(
			std::make_tuple(kHistList, monoHistList, infoList, lenList), size,
			threadNum);

	kHistList = nullptr;
	monoHistList = nullptr;
	infoList = nullptr;
	lenList = nullptr;

	isDataCleared = true;
}

template<class V>
MeanShift<V>::~MeanShift() {
	if (!isDataCleared) {
		clearData();
	}

	if (kHistList != nullptr || monoHistList != nullptr || infoList != nullptr
			|| lenList != nullptr) {
		cerr << "~MeanShift(): Data is not cleared!" << endl;
	}

	int s = clusterList->size();
#pragma omp parallel for schedule(static) num_threads(threadNum)
	for (int i = 0; i < s; i++) {
		delete clusterList->at(i);
	}
	delete clusterList;

	if (assignList != nullptr) {
		delete[] assignList;
	}
}

template<class V>
void MeanShift<V>::initData(Block *block) {
	size = block->size();
	auto tup = identity.unpackBlock(block, threadNum);
	kHistList = std::get < 0 > (tup);
	monoHistList = std::get < 1 > (tup);
	infoList = std::get < 2 > (tup);
	lenList = std::get < 3 > (tup);

	isDataCleared = false;
}

template<class V>
void MeanShift<V>::start() {

	//	if (threshold >= 0.9) {
	//		// double slackMax = Parameters::getMsSlackMax();
	//		//double slack = identity.getError() > slackMax ? slackMax : identity.getError();
	//		// threshold = threshold - slack;
	//		threshold = threshold - identity.getError();
	//		std::cout << "Relaxed threshold is: " << threshold << std::endl;
	//	}

	if (threshold <= 0.7) {
		isLowIdentity = true;
		//		double slackMax = Parameters::getMsSlackMax();
		//		double slack = identity.getError() > slackMax ? slackMax : identity.getError();
		//		// Relaxation is capped at 2% max
		//		threshold = threshold - slack;
	} else {
		isLowIdentity = false;
	}

	mergeThreshold = threshold;
	// Relaxation is not capped
	threshold = threshold - identity.getError();

	isMatrixCleared = false;
	initClusters();
	run(Parameters::getMsItr(), true);
	removeEmpty();
}
/**
 * !!! This method is the bottleneck for memory, especially on large genomes !!!
 */
template<class V>
void MeanShift<V>::initClusters() {
//cout << "Bottleneck: Before all-vs-all matrix." << endl;
	a = identity.score(kHistList, monoHistList, size, threadNum, lenList);
//cout << "Bottleneck: After all-vs-all matrix." << endl;
	// Start with every sequence as a cluster
	clusterList = new std::vector<Cluster<V>*>(size, nullptr);
	int kHistSize = identity.getKHistSize();
	int monoHistSize = identity.getMonoHistSize();

#pragma omp parallel for schedule(static) num_threads(threadNum)
	for (int i = 0; i < size; i++) {
		clusterList->at(i) = new Cluster<V>(kHistList, monoHistList, &a(i, 0),
				size, kHistSize, monoHistSize, threshold, i);
	}
//cout << "Bottleneck: After initializing all clusters." << endl;
}

template<class V>
void MeanShift<V>::run(int itrNum, bool canAssign) {
	if (itrNum > 1 && !canAssign) {
		cerr << "MeanShift<V>::run canAssign must be true ";
		cerr << "when the number of iterations is greater than 1." << endl;
		throw std::exception();
	}

	int count = 0;
	int oldClusterNumber = clusterList->size();
	//int itrNum = Parameters::getMsItr();
	int i;
	for (i = 0; i < itrNum; i++) {
		shift();
		if (isLowIdentity) {
			selectRepresentative();
			mergeGreedy(); // merge on representative sequences
		} else {
			mergeGreedy(); // merge on the synthetic means
			selectRepresentative();
		}

		if (canAssign) {
			updateIdentityList();
		}

		// Check convergence
		int newClusterNumber = clusterList->size();
		if (oldClusterNumber == newClusterNumber) {
			count++;
		} else {
			count = 0;
		}

		if (count == 2) {
			break;
		}
		oldClusterNumber = newClusterNumber;

		// The 'a' matrix should be deleted as soon as possible
		if (!isMatrixCleared) {
			copyIdentityList();
			a.clearData();
			isMatrixCleared = true;
		}
	}

	// Assign points to clusters -> contribution increases potentially
	if (canAssign) {
		assign();
	}

	// Once converged update members to include the new ones in the cluster
	// ToDo: Move this outside of run OR follow the same pattern as assign()
	// updateAccumulatedMean();

	// cout << "Converged in " << i << " iterations." << endl;
}

//template<class V>
//void MeanShift<V>::runOnce() {
//	if (isLowIdentity) {
//		shift();
//		selectRepresentative();
//		mergeGreedy(); // merge on representative sequences
//	} else {
//		shift();
//		mergeGreedy(); // merge on the synthetic means
//		selectRepresentative();
//	}
//	// Once converged update members to include the new ones in the cluster
//	updateAccumulatedMean();
//}

template<class V>
void MeanShift<V>::shift() {
	int s = clusterList->size();
#pragma omp parallel for schedule(static) num_threads(threadNum)
	for (int i = 0; i < s; i++) {
		clusterList->at(i)->shiftWeighted();
	}
}

template<class V>
void MeanShift<V>::selectRepresentative() {
	int s = clusterList->size();
#pragma omp parallel for schedule(static) num_threads(threadNum)
	for (int i = 0; i < s; i++) {
		auto c = clusterList->at(i);
		if (c->getHasShifted()) {

			// Calculate identities between members and the new center
			auto l = c->getMemberList();
			int m = l->size();
			V *kList[m];
			uint64_t *mList[m];
			int lList[m];
			for (int x = 0; x < m; x++) {
				int z = l->at(x);
				kList[x] = kHistList[z];
				mList[x] = monoHistList[z];
				lList[x] = lenList[z];
			}

			double *v = identity.score(c->getKHistMean(), kList,
					c->getMonoHistMean(), mList, m, 1 /*threadNum*/,
					c->getLength(), lList);

			//c->clearMemberList();

			// Calculate identity between the old center and the new
			double h = -1;
			if (c->getKHistOld() != nullptr && c->getMonoHistOld() != nullptr) {
				h = identity.score(c->getKHistMean(), c->getKHistOld(),
						c->getMonoHistMean(), c->getMonoHistOld(),
						identity.calcRatio(c->getLength(), c->getOldLength()),
						c->getLength(), c->getOldLength());
			}

			// Find the representative with the highest identity with the new center
			double max = h;
			int index = -1;
			for (int e = 0; e < m; e++) {
				if (v[e] > max) {
					max = v[e];
					index = e;
				}
			}

			// Set a new rep or keep the old one
			if (index >= 0) {
				c->setRepresentative(kList[index], mList[index], false);
			} else {
				c->setRepresentative(c->getKHistOld(), c->getMonoHistOld(),
						true);
			}

			// Clean up
			delete[] v;
		}
	}
}

/**
 * This method finds the center closest to the mean of all centers to
 * be merged and returns it. Also, this center is removed from the list
 *
 * Currently unused
 *
 */
//template<class V>
//Cluster<V>* MeanShift<V>::mergeGreedyHelper(std::vector<Cluster<V>*> &list) {
//	Cluster < V > *clusterToMerge = nullptr;
//	if (list.empty()) {
//		cerr << "mergeGreedyHelper is expecting a non-empty vector." << endl;
//		throw std::exception();
//	} else if (list.size() == 1) {
//		clusterToMerge = list[0];
//	} else {
//		double size = list.size();
//
//		V *kList[(int) size];
//		uint64_t *monoList[(int) size];
//		int lList[(int) size];
//
//		int kHistSize = identity.getKHistSize();
//		int monoHistSize = identity.getMonoHistSize();
//
//		V *kMean = new V[kHistSize] { 0 };
//		uint64_t monoMean[monoHistSize] { 0 };
//
//		for (int i = 0; i < size; i++) {
//			auto kHist = list[i]->getKHistMean();
//			auto monoHist = list[i]->getMonoHistMean();
//
//			kList[i] = kHist;
//			monoList[i] = monoHist;
//			lList[i] = list[i]->getLength();
//
//			for (int j = 0; j < kHistSize; j++) {
//				kMean[j] += kHist[j];
//			}
//
//			for (int j = 0; j < monoHistSize; j++) {
//				monoMean[j] += monoHist[j];
//			}
//		}
//
//		for (int j = 0; j < kHistSize; j++) {
//			kMean[j] = round(kMean[j] / size);
//		}
//
//		int len = 0;
//		for (int j = 0; j < monoHistSize; j++) {
//			monoMean[j] = round(monoMean[j] / size);
//			len += monoMean[j];
//		}
//
//		// Calculate the identity scores using identity
//		double *idList = identity.score(kMean, kList, monoMean, monoList, size,
//				threadNum, len, lList);
//
//		int index = std::max_element(idList, idList + (int) size) - idList;
//		if (index < 0 || index >= size) {
//			cerr << "mergeGreedyHelper produced invalid index of " << index
//					<< " while the size is " << (int) size << endl;
//		}
//		cout << "Max index is: " << index << endl;
//
//		clusterToMerge = list.at(index);
//		list.erase(list.begin() + index);
//
//		delete[] kMean;
//	}
//	return clusterToMerge;
//}

/**
 * The merge step is first come first served
 *  1. Start with the first cluster
 * 	2. Find all clusters that are close to it
 * 	3. Merge the largest cluster among these with the rest
 * 	4. Remove and Delete these clusters
 *  5. Continue with a new cluster
 */
template<class V>
void MeanShift<V>::mergeGreedy() {
	int s = clusterList->size();
	std::vector<Cluster<V>*> *mergedList = new std::vector<Cluster<V>*>();

	// If an entry is true, the cluster is not merged
	std::vector<bool> remainList(s, true);

	for (int i = 0; i < s; i++) {
		// Find a cluster to merge
		if (remainList[i]) {
			// Mark this cluster done
			remainList[i] = false;
			auto cluster = clusterList->at(i);
			mergedList->push_back(cluster);

			// Collect indexes of remaining clusters
			std::vector<int> indexList;
			for (int h = i + 1; h < s; h++) {
				if (remainList[h]) {
					indexList.push_back(h);
				}
			}
			int r = indexList.size();

			if (r > 0) {
				// Collect means (kmer and monomer) in two arrays
				V *kList[r];
				uint64_t *monoList[r];
				int cLenList[r];

#pragma omp parallel for schedule(static) num_threads(threadNum)
				for (int y = 0; y < r; y++) {
					auto c = clusterList->at(indexList[y]);
					kList[y] = c->getKHistMean();
					monoList[y] = c->getMonoHistMean();
					cLenList[y] = c->getLength();
				}

				double *idList = identity.score(cluster->getKHistMean(), kList,
						cluster->getMonoHistMean(), monoList, r, threadNum,
						cluster->getLength(), cLenList);

				std::vector<Cluster<V>*> similarList;
				for (int u = 0; u < r; u++) {
					//if (idList[u] >= threshold) {
					if (idList[u] >= mergeThreshold) {
						similarList.push_back(clusterList->at(indexList[u]));
						remainList[indexList[u]] = false;
					}
				}
				if (similarList.size() > 0) {
					// cluster->mergeWeighted(similarList);
					// cluster->mergeLargest(similarList);

					// similarList.push_back(cluster);
					// mergedList->back() = mergeGreedyHelper(similarList);
					// mergedList->back()->mergeSimple(similarList);

					similarList.push_back(cluster);
					int max = similarList[0]->getContribution();
					int index = 0;
					int simCount = similarList.size();
					for (int a = 1; a < simCount; a++) {
						int w = similarList[a]->getContribution();
						if (w > max) {
							max = w;
							index = a;
						}
					}
					mergedList->back() = similarList[index];
					similarList.erase(similarList.begin() + index);
					mergedList->back()->mergeSimple(similarList);

					for (auto d : similarList) {
						delete d;
					}
				}
				delete[] idList;
			} else {
				break;
			}
		}
	}

	if (clusterList->size() < mergedList->size()) {
		std::cerr << "Merged list cannot be larger than original list.";
		std::cerr << std::endl;
		std::cerr << mergedList->size() << " vs. " << clusterList->size();
		std::cerr << std::endl;
		throw std::exception();
	}

	delete clusterList;
	clusterList = mergedList;
}

/**
 * Merge each connected component
 */
//template<class V>
//void MeanShift<V>::merge() {
//	int s = clusterList->size();
//	Matrix ava = calcAllCenterVsAllCenter();
//
//	auto p = ClusteringUtil::findConnectedComponents(ava, mergeThreshold);
//	// auto p = ClusteringUtil::findConnectedComponents(ava, threshold);
//	vector<int> flagList = p.first;
//	int compNum = p.second;
//
//// Separate clusters of each component in a vector
//	vector < vector<Cluster<V>*> > compList;
//	compList.reserve(compNum);
//	for (int i = 0; i < compNum; i++) {
//		compList.push_back(vector<Cluster<V>*>());
//	}
//	for (int i = 0; i < s; i++) {
//		compList[flagList[i] - 1].push_back(clusterList->at(i));
//	}
//
//// Clear old clusters
//	clusterList->clear();
//
//// Merge clusters in each component
//	for (int i = 0; i < compNum; i++) {
//		// Find the cluster with the largest w
//		double w = -1;
//		int index = -1;
//		int size = compList[i].size();
//		for (int j = 0; j < size; j++) {
//			Cluster < V > *c = compList[i][j];
//			if (c->getContribution() > w) {
//				w = c->getContribution();
//				index = j;
//			}
//		}
//
//		if (index == -1) {
//			std::cerr << "Could not find the cluster with the largest w";
//			std::cerr << std::endl;
//			throw std::exception();
//		}
//
//		// Put the largest cluster at the back
//		if (index != size - 1) {
//			Cluster < V > *temp = compList[i][size - 1];
//			compList[i][size - 1] = compList[i][index];
//			compList[i][index] = temp;
//			//std::cout << compList[i].back()->getContribution() << std::endl;
//		}
//
//		auto c = compList[i].back();
//		compList[i].pop_back();
//		// c->mergeNew(compList[i]); // ToDo
//		c->mergeWeighted(compList[i]);
//		clusterList->push_back(c);
//		for (auto x : compList[i]) {
//			delete x;
//		}
//	}
//}

/**
 * This method should be called on the clusters remaining when
 * there is a new center, i.e., a new representative sequence
 */
template<class V>
void MeanShift<V>::updateIdentityList() {
	int n = clusterList->size();
	for (int i = 0; i < n; i++) {
		Cluster < V > *c = clusterList->at(i);

		// Calculate identities between a block and the new center when
		// there is a new representative
		if (!c->getIsIdListUpToDate()) {
			double *v = identity.score(c->getKHistMean(), kHistList,
					c->getMonoHistMean(), monoHistList, size, threadNum,
					c->getLength(), lenList);

			c->setIdentityList(v);
		}
	}
}

template<class V>
void MeanShift<V>::copyIdentityList() {
	int n = clusterList->size();
#pragma omp parallel for schedule(static) num_threads(threadNum)
	for (int i = 0; i < n; i++) {
		auto cluster = clusterList->at(i);
		if (cluster->getIsOriginalIdList()) {
			cluster->copyIdentityList();
		}
	}
}

template<class V>
void MeanShift<V>::updateAccumulatedMean() {
	int n = clusterList->size();
#pragma omp parallel for schedule(static) num_threads(threadNum)
	for (int i = 0; i < n; i++) {
		clusterList->at(i)->updateAccumulatedMean();
	}
}

template<class V>
void MeanShift<V>::assign() {
	if (assignList != nullptr) {
		delete[] assignList;
	}

	// Which cluster a point belongs to
	assignList = new int[size];
	std::fill_n(assignList, size, -1);

	// Best identity score
	vector<double> scoreList(size, -1.0);

	int clusterNum = clusterList->size();
	for (int i = 0; i < clusterNum; i++) {
		auto cluster = clusterList->at(i);
		const double *identityList = cluster->getIdentityList();
#pragma omp parallel for schedule(static) num_threads(threadNum)
		for (int j = 0; j < size; j++) {
			double score = identityList[j];
			if (score >= threshold && score > scoreList[j]) {
				assignList[j] = i;
				scoreList[j] = score;
			}
		}
	}

	for (int i = 0; i < size; i++) {
		if (assignList[i] > -1) {
			clusterList->at(assignList[i])->incrementAssignment();
		}
	}
}

/**
 * Running the mean shift using this method does not result in new clusters
 */
template<class V>
void MeanShift<V>::updateReferenceData(Block *block) {
	//clearData();
	initData(block);
	int clusterNum = clusterList->size();
#pragma omp parallel for schedule(static) num_threads(threadNum)
	for (int i = 0; i < clusterNum; i++) {
		clusterList->at(i)->updateReferenceData(kHistList, monoHistList, size);
	}

	updateIdentityList();
}

/**
 * Add additional clusters from previous or concurrent runs
 */
template<class V>
void MeanShift<V>::addClusters(const vector<Cluster<V>*> *otherClusterList) {
	int kHistSize = identity.getKHistSize();
	int monoHistSize = identity.getMonoHistSize();

	for (auto c : *otherClusterList) {
		auto kHist = c->getKHistMean();
		auto monoHist = c->getMonoHistMean();
		double *idList = identity.score(kHist, kHistList, monoHist,
				monoHistList, size, threadNum, c->getLength(), lenList);

		auto h = new Cluster<V>(kHistList, monoHistList, idList, size,
				kHistSize, monoHistSize, threshold, kHist, monoHist,
				c->getContribution(), c->getAssignment());
		clusterList->push_back(h);
	}
}

template<class V>
Matrix MeanShift<V>::calcAllCenterVsAllCenter() {
	int s = clusterList->size();

	if (s > Parameters::getMsMaxMatrixSize()) {
		std::cerr << "MeanShift: Matrix is too large." << std::endl;
		throw std::exception();
	}

	// Collect means (kmer and monomer) in two arrays
	V *kList[s];
	uint64_t *monoList[s];
	int cLenList[s];
#pragma omp parallel for schedule(static) num_threads(threadNum)
	for (int i = 0; i < s; i++) {
		auto c = clusterList->at(i);
		kList[i] = c->getKHistMean();
		monoList[i] = c->getMonoHistMean();
		cLenList[i] = c->getLength();
	}

	// Perform all versus all
	Matrix ava = identity.score(kList, monoList, s, threadNum, cLenList);

	return ava;
}

/**
 * Remove clusters that did not merge at all, i.e. single sequences
 */
template<class V>
void MeanShift<V>::removeSingles() {
	vector<Cluster<V>*> *newList = new vector<Cluster<V>*>();
	newList->reserve(clusterList->size());
	for (auto c : *clusterList) {
		if (c->getContribution() > 1) {
			newList->push_back(c);
		} else {
			delete c;
		}
	}

	delete clusterList;
	clusterList = newList;
}

/**
 * Remove cluster with no points assigned to
 */
template<class V>
void MeanShift<V>::removeEmpty() {
	vector<Cluster<V>*> *newList = new vector<Cluster<V>*>();
	newList->reserve(clusterList->size());
	for (auto c : *clusterList) {
		if (c->getAssignment() > 1) {
			newList->push_back(c);
		} else {
			delete c;
		}
	}

	delete clusterList;
	clusterList = newList;
}
/**
 * This method return copies of unassigned data points.
 * It is the client responsibility to free the memory allocated to these points.
 *
 */
template<class V>
std::tuple<V**, uint64_t**, std::string**, int*, int> MeanShift<V>::findUnassignedData() {
	if (assignList == nullptr) {
		std::cerr << "The assignment list is null. ";
		std::cerr << "Please run the mean shift first." << std::endl;
		throw std::exception();
	}

	// Count unassigned data points
	vector<int> indexList;
	indexList.reserve(size);
	for (int j = 0; j < size; j++) {
		if (assignList[j] == -1) {
			indexList.push_back(j);
		}
	}
	int unassignedCount = indexList.size();

	// Copy unassigned data points
	int kHistSize = identity.getKHistSize();
	int monoHistSize = identity.getMonoHistSize();

	V **kHistListU = new V*[unassignedCount];
	uint64_t **monoHistListU = new uint64_t*[unassignedCount];
	std::string **infoListU = new std::string*[unassignedCount];
	int *lenListU = new int[unassignedCount];

#pragma omp parallel for schedule(static) num_threads(threadNum)
	for (int i = 0; i < unassignedCount; i++) {
		int j = indexList[i];

		// ! For the std::copy the end is exclusive
		V *kHist = new V[kHistSize];
		std::copy(&kHistList[j][0], &kHistList[j][kHistSize], kHist);
		kHistListU[i] = kHist;

		uint64_t *monoHist = new uint64_t[monoHistSize];
		std::copy(&monoHistList[j][0], &monoHistList[j][monoHistSize],
				monoHist);
		monoHistListU[i] = monoHist;

		infoListU[i] = new string(*infoList[j]);

		lenListU[i] = lenList[j];
	}

	return std::make_tuple(kHistListU, monoHistListU, infoListU, lenListU,
			unassignedCount);
}

template<class V>
const vector<Cluster<V>*>* MeanShift<V>::getClusterList() const {
	return clusterList;
}

//template<class V>
//double MeanShift<V>::getThreshold() {
//	return threshold;
//}

//template<class V>
//void MeanShift<V>::checkClusters() {
//	int num = clusterList->size();
//	for (int i = 0; i < num; i++) {
//		if (clusterList->at(i)->isAllZeros()) {
//			std::cerr << "Found a center made of all zeros." << std::endl;
//			throw std::exception();
//		}
//	}
//}
