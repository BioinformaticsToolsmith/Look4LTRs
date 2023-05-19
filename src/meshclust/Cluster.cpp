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
 * Cluster.cpp
 *
 *  Created on: Dec 23, 2020
 *      Author: Hani Z. Girgis, PhD
 */

/**
 * kList: A list of k-mer histograms
 * monoList: A list of monomer histograms
 * idList: A list of identity scores of the initial center to all other sequences
 * Memory: This class is responsible for freeing memory used by idList
 * lSize: Number of sequences in a set
 * kSize: The size of a k-mer histogram
 * monoSize: The size of a monomer histogram
 * index: is the index of the initial center of this cluster
 *
 * !!!!Sometime the identity list is a row in the identity matrix.
 * If this is the case it will not be deleted.
 */
template<class V>
Cluster<V>::Cluster(V **kList, uint64_t **monoList, double *idList, int lSize,
		int kSize, int monoSize, double t, int index) {

	init(kList, monoList, idList, lSize, kSize, monoSize, t);

	memcpy(kHistMean, kList[index], kHistSize * sizeof *kHistMean);
	// memcpy(kHistOld, kList[index], kHistSize * sizeof *kHistMean);
	memcpy(monoHistMean, monoList[index], monoHistSize * sizeof *monoHistMean);
	// memcpy(monoHistOld, monoList[index], monoHistSize * sizeof *monoHistMean);

	//canDelete = false;
	isOriginalIdList = true;
	contribution = 0;
	oldContribution = 0;
	assignment = 0;
}

template<class V>
Cluster<V>::Cluster(V **kList, uint64_t **monoList, double *idList, int lSize,
		int kSize, int monoSize, double t, V *kHist, uint64_t *monoHist,
		int count, int assign) {

	init(kList, monoList, idList, lSize, kSize, monoSize, t);

	kHistOld = new V[kHistSize] { 0 };
	monoHistOld = new uint64_t[monoHistSize] { 0 };

	memcpy(kHistMean, kHist, kHistSize * sizeof *kHistMean);
	memcpy(kHistOld, kHist, kHistSize * sizeof *kHistMean);
	memcpy(monoHistMean, monoHist, monoHistSize * sizeof *monoHistMean);
	memcpy(monoHistOld, monoHist, monoHistSize * sizeof *monoHistMean);

	contribution = count;
	oldContribution = count;
	assignment = assign;
	isOriginalIdList = false;
}

template<class V>
void Cluster<V>::init(V **kList, uint64_t **monoList, double *idList, int lSize,
		int kSize, int monoSize, double t) {
	kHistList = kList;
	monoHistList = monoList;
	listSize = lSize;
	kHistSize = kSize;
	monoHistSize = monoSize;
	identityList = idList;
	threshold = t;
	isIdListUpToDate = true;
	hasShifted = true;

	// Elements of these two arrays are initialized in the first call
	// to shift if they are not provided
	kHistMean = new V[kHistSize] { 0 };
	monoHistMean = new uint64_t[monoHistSize] { 0 };

	kHistOld = nullptr;
	monoHistOld = nullptr;

	memberList = new std::vector<int>();
}

template<class V>
Cluster<V>::~Cluster() {
	if (!isOriginalIdList) {
		delete[] identityList;
	}
	delete[] kHistMean;
	delete[] monoHistMean;

	if (kHistOld != nullptr) {
		delete[] kHistOld;
	}
	if (monoHistOld != nullptr) {
		delete[] monoHistOld;
	}

	delete memberList;
}

/**
 * Equation 8 in Cheng 1995
 *
 * Notes:
 * 1. This method assumes that the identity scores between the center of
 * the cluster and every other point is up to date: Enforced
 * 2. This method applies flat kernel: x >= threshold? 1.0 : 0.0
 */
template<class V>
void Cluster<V>::shiftWeighted() {
	// Make sure that the identity list is up to date
	if (!isIdListUpToDate) {
		cerr << "shiftWeighted is expecting an up-to-date identity list.";
		cerr << endl;
		throw std::exception();
	}

	double *kTemp = new double[kHistSize] { 0.0 }; // May be a large array
	double monoTemp[monoHistSize] { 0.0 };

	double n = 0.0; // Number of new sequences within the threshold
	memberList->clear();

	for (int i = 0; i < listSize; i++) {
		if (identityList[i] >= threshold) {
			n++;
			memberList->push_back(i);

			V *kSeqHist = kHistList[i];
			uint64_t *monoSeqHist = monoHistList[i];

			for (int j = 0; j < kHistSize; j++) {
				kTemp[j] += kSeqHist[j];
			}

			for (int j = 0; j < monoHistSize; j++) {
				monoTemp[j] += monoSeqHist[j];
			}
		}
	}

	// If no points nearby this center, it does not shift
	if (n >= 1.0) {
		if (oldContribution > 1) {
			if (kHistOld == nullptr || monoHistOld == nullptr) {
				cerr << "Cluster<V>::shiftWeighted(): "
						<< "one of the old histograms is null.";
				cerr << endl;
				throw std::exception();
			}

			double total = oldContribution + n;
			double oW = oldContribution / total;
			double nW = n / total;
			double nWByN = nW / n;

			for (int j = 0; j < kHistSize; j++) {
				kHistMean[j] = round((nWByN * kTemp[j]) + (oW * kHistOld[j]));
			}

			for (int j = 0; j < monoHistSize; j++) {
				monoHistMean[j] = round(
						(nWByN * monoTemp[j]) + (oW * monoHistOld[j]));
			}
		} else {
			for (int j = 0; j < kHistSize; j++) {
				kHistMean[j] = round(kTemp[j] / n);
			}

			for (int j = 0; j < monoHistSize; j++) {
				monoHistMean[j] = round(monoTemp[j] / n);
			}
		}
		hasShifted = true;
	} else {
		hasShifted = false;
	}

	contribution = n + oldContribution;

	delete[] kTemp;
}

/**
 * The shifted
 */
template<class V>
void Cluster<V>::setRepresentative(V *kHist, uint64_t *monoHist,
		bool isUpToDate) {
//	for (int i = 0; i < kHistSize; i++) {
//		kHistMean[i] = kHist[i];
//
//	}
//	for (int j = 0; j < monoHistSize; j++) {
//		monoHistMean[j] = monoHist[j];
//	}

	memcpy(kHistMean, kHist, kHistSize * sizeof *kHistMean);
	memcpy(monoHistMean, monoHist, monoHistSize * sizeof *monoHistMean);

	isIdListUpToDate = isUpToDate;
}

template<class V>
void Cluster<V>::mergeSimple(std::vector<Cluster*> &clusterList) {
	if (!clusterList.empty()) {
		// Perform union on members contributed to all merged clusters
		std::unordered_set<int> s;
		s.insert(memberList->begin(), memberList->end());

		for (auto c : clusterList) {
			contribution += c->getContribution();
			oldContribution += c->getOldContribution();
			assignment += c->getAssignment();
			auto l = c->getMemberList();
			s.insert(l->begin(), l->end());
		}

		memberList->clear();
		memberList->reserve(s.size());
		memberList->insert(memberList->begin(), s.begin(), s.end());

//#pragma omp critical
//		{
//			for (auto c : clusterList) {
//				cout << c->getMemberList()->size() << endl;
//			}
//			cout << "Set size: " << s.size() << " list size: "
//					<< memberList->size() << endl;
//			cout << "-----" << endl;
//		}
	}
}

template<class V>
vector<int>* Cluster<V>::getMemberList() const {
	return memberList;
}

template<class V>
void Cluster<V>::clearMemberList() {
	memberList->clear();
}

/**
 * Currently unused
 */
//template<class V>
//void Cluster<V>::mergeWeightedHelper(std::vector<Cluster*> &clusterList,
//		V *kMean, uint64_t *monoMean) {
//	// Sum all contributions
//	// Add pseudo count to handle the initial clusters
//	double total = contribution + 1;
//	for (auto c : clusterList) {
//		total += (c->getContributionCount() + 1);
//	}
//
//	// Initialize the new means according to my weight
//	double *kTemp = new double[kHistSize]; // May be a large array
//	double monoTemp[monoHistSize];
//
//	double myWeight = (contribution + 1) / total;
//	for (int j = 0; j < kHistSize; j++) {
//		kTemp[j] = myWeight * kMean[j];
//	}
//
//	for (int j = 0; j < monoHistSize; j++) {
//		monoTemp[j] = myWeight * monoMean[j];
//	}
//
//	// Add each mean according to its weight
//	for (auto c : clusterList) {
//		double otherWeight = (c->getContributionCount() + 1) / total;
//		V *kHist = c->getKHistMean();
//		uint64_t *monoHist = c->getMonoHistMean();
//
//		for (int j = 0; j < kHistSize; j++) {
//			kTemp[j] += otherWeight * kHist[j];
//		}
//
//		for (int j = 0; j < monoHistSize; j++) {
//			monoTemp[j] += otherWeight * monoHist[j];
//		}
//	}
//
//	// Round every entry
//	for (int j = 0; j < kHistSize; j++) {
//		kMean[j] = round(kTemp[j]);
//	}
//
//	for (int j = 0; j < monoHistSize; j++) {
//		monoMean[j] = round(monoTemp[j]);
//	}
//
//	delete[] kTemp;
//}
/**
 *
 * Currently unused
 *
 */
//template<class V>
//void Cluster<V>::mergeWeighted(std::vector<Cluster*> &clusterList) {
//	if (!clusterList.empty()) {
//		mergeWeightedHelper(clusterList, kHistMean, monoHistMean);
//
//		// Check
//		//		if (isAllZeros()) {
//		//			std::cerr << "Found a center made of all zeros." << std::endl;
//		//			throw std::exception();
//		//		}
//
//		mergeWeightedHelper(clusterList, kHistOld, monoHistOld);
//
//		// Update member count
//		for (auto c : clusterList) {
//			contribution += c->getContribution();
//			oldContribution += c->getOldContribution();
//			assignment += c->getAssignment();
//		}
//	}
//}
/**
 *
 * Currently unused
 *
 */
//template<class V>
//void Cluster<V>::mergeLargest(std::vector<Cluster*> &clusterList) {
//	if (!clusterList.empty()) {
//		//mergeWeightedHelper(clusterList, kHistMean, monoHistMean);
//		//mergeWeightedHelper(clusterList, kHistOld, monoHistOld);
//
//		int clustNum = clusterList.size();
//		int index = -1;
//		int max = getContribution();
//		for (int u = 0; u < clustNum; u++) {
//			int c = clusterList[u]->getContribution();
//			if (c > max) {
//				max = c;
//				index = u;
//			}
//		}
//
//		if (index != -1) {
//			auto o = clusterList[index];
//			V *oKHist = o->getKHistMean();
//			uint64_t *oMonoHist = o->getMonoHistMean();
//
//			V *oKHistOld = o->getKHistOld();
//			uint64_t *oMonoHistOld = o->getMonoHistOld();
//
//			std::copy(oKHist, oKHist + kHistSize, kHistMean);
//			std::copy(oMonoHist, oMonoHist + monoHistSize, monoHistMean);
//			if (oKHistOld != nullptr && kHistOld != nullptr) {
//				std::copy(oKHistOld, oKHistOld + kHistSize, kHistOld);
//			}
//			if (oMonoHistOld != nullptr && monoHistOld != nullptr) {
//				std::copy(oMonoHistOld, oMonoHistOld + monoHistSize,
//						monoHistOld);
//			}
//		}
//
//		// Update member count
//		for (auto c : clusterList) {
//			contribution += c->getContribution();
//			oldContribution += c->getOldContribution();
//			assignment += c->getAssignment();
//		}
//	}
//}

/**
 * Note: This method frees the memory allocated to the old identity list
 */
template<class V>
void Cluster<V>::setIdentityList(double *newList) {
	if (!isOriginalIdList) {
		delete[] identityList;
	}

	identityList = newList;
	isIdListUpToDate = true;
	isOriginalIdList = false;
}

template<class V>
void Cluster<V>::copyIdentityList() {
	if (!isOriginalIdList) {
		cerr << "copyIdentityList can be applied only to ";
		cerr << "a row of the all-vs-all matrix." << endl;
		throw std::exception();
	}

	double *copyList = new double[listSize];
	memcpy(copyList, identityList, listSize * sizeof *identityList);

	identityList = copyList;
	isOriginalIdList = false;
}

template<class V>
void Cluster<V>::updateAccumulatedMean() {
	// Update old mean, i.e., mean constructed from all previous blocks
	if (kHistOld == nullptr) {
		kHistOld = new V[kHistSize];
	}
	if (monoHistOld == nullptr) {
		monoHistOld = new uint64_t[monoHistSize];
	}

	memcpy(kHistOld, kHistMean, kHistSize * sizeof *kHistMean);
	memcpy(monoHistOld, monoHistMean, monoHistSize * sizeof *monoHistMean);

	oldContribution = contribution;
}

template<class V>
void Cluster<V>::updateReferenceData(V **kList, uint64_t **monoList,
		int lSize) {
	kHistList = kList;
	monoHistList = monoList;
	listSize = lSize;
	isIdListUpToDate = false;
}

// OK
template<class V>
V* Cluster<V>::getKHistMean() const {
	return kHistMean;
}

// OK
template<class V>
uint64_t* Cluster<V>::getMonoHistMean() const {
	return monoHistMean;
}

// OK
template<class V>
V* Cluster<V>::getKHistOld() const {
	return kHistOld;
}

// OK
template<class V>
uint64_t* Cluster<V>::getMonoHistOld() const {
	return monoHistOld;
}

template<class V>
const double* Cluster<V>::getIdentityList() {
	return identityList;
}

template<class V>
int Cluster<V>::getContribution() const {
	return contribution;
}

template<class V>
int Cluster<V>::getOldContribution() const {
	return oldContribution;
}

// OK
template<class V>
int Cluster<V>::getLength() const {
	return std::accumulate(monoHistMean, monoHistMean + monoHistSize, 0,
			std::plus<int>());
}

// OK
template<class V>
int Cluster<V>::getOldLength() const {
	return std::accumulate(monoHistOld, monoHistOld + monoHistSize, 0,
			std::plus<int>());
}

template<class V>
void Cluster<V>::incrementAssignment() {
	assignment++;
}

template<class V>
int Cluster<V>::getAssignment() {
	return assignment;
}

template<class V>
bool Cluster<V>::getIsIdListUpToDate() const {
	return isIdListUpToDate;
}

template<class V>
bool Cluster<V>::getHasShifted() const {
	return hasShifted;
}

template<class V>
bool Cluster<V>::getIsOriginalIdList() const {
	return isOriginalIdList;
}

// Currently not used, but keep it
//template<class V>
//bool Cluster<V>::isAllZeros() {
//	bool r = Util::isAllZeros < V > (kHistMean, kHistSize)
//			|| Util::isAllZeros < uint64_t > (monoHistMean, monoHistSize);
//	return r;
//}
