cmake_minimum_required (VERSION 3.10)
project (identity)

# Credit: https://stackoverflow.com/questions/52180281/cmake-cxx-compiler-version-is-pointing-to-the-old-gcc-version
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 7.5)
        message(FATAL_ERROR "GNU g++ 7.5.0 or later is required. Your current version is: " ${CMAKE_CXX_COMPILER_VERSION})
    endif()  
else()
    message(FATAL_ERROR "Your compiler is currently unsupported: " ${CMAKE_CXX_COMPILER_ID})
endif()

set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})
set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR})
set(CMAKE_CXX_STANDARD_REQUIRED on)
set(CMAKE_CXX_STANDARD 17)
# -ggdb
# -O3 -ffast-math
# -g -O0
# -O2 -fgcse-after-reload -fipa-cp-clone -floop-interchange -floop-unroll-and-jam -fpeel-loops -fpredictive-commoning -fsplit-loops -fsplit-paths -ftree-loop-distribution -ftree-partial-pre -funswitch-loops -fvect-cost-model=dynamic -fversion-loops-for-strides -ffast-math
# There is a bug if -O3 if used in GCC 11.1.0, so we instead added -O2 with all options from -O3
set(CMAKE_CXX_FLAGS " -O2 -fgcse-after-reload -fipa-cp-clone -floop-interchange -floop-unroll-and-jam -fpeel-loops -fpredictive-commoning -fsplit-loops -fsplit-paths -ftree-loop-distribution -ftree-partial-pre -funswitch-loops -fvect-cost-model=dynamic -fversion-loops-for-strides -ffast-math  -fopenmp -fmessage-length=0 -fpermissive")

add_library(main SHARED		
${CMAKE_SOURCE_DIR}/src/FastaReader.cpp		
${CMAKE_SOURCE_DIR}/src/Mutator.cpp			
${CMAKE_SOURCE_DIR}/src/ReaderAlignerCoordinator.cpp			
${CMAKE_SOURCE_DIR}/src/Parameters.cpp			
${CMAKE_SOURCE_DIR}/src/Matrix.cpp
${CMAKE_SOURCE_DIR}/src/DataGenerator.cpp
${CMAKE_SOURCE_DIR}/src/SynDataGenerator.cpp
${CMAKE_SOURCE_DIR}/src/GLMClassifier.cpp
${CMAKE_SOURCE_DIR}/src/GLMRegressor.cpp
${CMAKE_SOURCE_DIR}/src/Normalizer.cpp
${CMAKE_SOURCE_DIR}/src/SimConverter.cpp
${CMAKE_SOURCE_DIR}/src/FeatureExpander.cpp
${CMAKE_SOURCE_DIR}/src/GLM.cpp
${CMAKE_SOURCE_DIR}/src/Node.cpp
${CMAKE_SOURCE_DIR}/src/Evaluator.cpp
${CMAKE_SOURCE_DIR}/src/Feature.cpp
${CMAKE_SOURCE_DIR}/src/FeatureSquared.cpp
${CMAKE_SOURCE_DIR}/src/FeaturePaired.cpp
${CMAKE_SOURCE_DIR}/src/StatisticInfo.cpp
${CMAKE_SOURCE_DIR}/src/ITransformer.cpp
${CMAKE_SOURCE_DIR}/src/GLMPredictor.cpp
${CMAKE_SOURCE_DIR}/src/Serializer.cpp
${CMAKE_SOURCE_DIR}/src/GlobAlignE.cpp
)

add_library(clustering SHARED
${CMAKE_SOURCE_DIR}/src/meshclust/MeShClust.cpp
${CMAKE_SOURCE_DIR}/src/meshclust/ClusterInfo.cpp
${CMAKE_SOURCE_DIR}/src/meshclust/ClusterEvaluator.cpp
${CMAKE_SOURCE_DIR}/src/meshclust/GMM.cpp
)
target_link_libraries(clustering main)

set(LTR_HEADER_FILES 
${CMAKE_SOURCE_DIR}/src/ltr/LtrParameters.h
${CMAKE_SOURCE_DIR}/src/ltr/Node.h
${CMAKE_SOURCE_DIR}/src/ltr/DirectedGraph.h
${CMAKE_SOURCE_DIR}/src/ltr/RT.h
${CMAKE_SOURCE_DIR}/src/ltr/DeepNesting.fwd.h
${CMAKE_SOURCE_DIR}/src/ltr/ModulePipeline.fwd.h
)

add_library(ltr SHARED
${CMAKE_SOURCE_DIR}/src/ltr/LtrParameters.cpp
${CMAKE_SOURCE_DIR}/src/ltr/Stretch.cpp
${CMAKE_SOURCE_DIR}/src/ltr/Merger.cpp
${CMAKE_SOURCE_DIR}/src/ltr/ScorerTr.cpp
${CMAKE_SOURCE_DIR}/src/ltr/StretchFeature.cpp
${CMAKE_SOURCE_DIR}/src/ltr/Element.cpp
${CMAKE_SOURCE_DIR}/src/ltr/RTComplete.cpp
${CMAKE_SOURCE_DIR}/src/ltr/RTSolo.cpp
${CMAKE_SOURCE_DIR}/src/ltr/RTComplex.cpp
${CMAKE_SOURCE_DIR}/src/ltr/Detector.cpp
${CMAKE_SOURCE_DIR}/src/ltr/SGD.cpp
${CMAKE_SOURCE_DIR}/src/ltr/Matcher.cpp
${CMAKE_SOURCE_DIR}/src/ltr/CaseMatcher.cpp
${CMAKE_SOURCE_DIR}/src/ltr/CaseSingle.cpp
${CMAKE_SOURCE_DIR}/src/ltr/CaseSolo.cpp
${CMAKE_SOURCE_DIR}/src/ltr/CaseRecentComplex.cpp
${CMAKE_SOURCE_DIR}/src/ltr/CaseRecent.cpp
${CMAKE_SOURCE_DIR}/src/ltr/Output.cpp
${CMAKE_SOURCE_DIR}/src/ltr/OutputBed.cpp
${CMAKE_SOURCE_DIR}/src/ltr/OutputRtr.cpp
${CMAKE_SOURCE_DIR}/src/ltr/OutputCpx.cpp
${CMAKE_SOURCE_DIR}/src/ltr/Filter.cpp
${CMAKE_SOURCE_DIR}/src/ltr/LtrUtility.cpp
${CMAKE_SOURCE_DIR}/src/ltr/LocalAlignment.cpp
${CMAKE_SOURCE_DIR}/src/ltr/PostProcess.cpp
${CMAKE_SOURCE_DIR}/src/ltr/ModulePipeline.cpp
${CMAKE_SOURCE_DIR}/src/ltr/DeepNesting.cpp
)
target_link_libraries(ltr main)

set(MSL_SOURCE
${CMAKE_SOURCE_DIR}/src/msl/MslToMatrix.cpp
)

set(TEST_SOURCE
${CMAKE_SOURCE_DIR}/src/test/Test.cpp
)

set(HEADER_FILES 
${CMAKE_SOURCE_DIR}/src/Util.h
${CMAKE_SOURCE_DIR}/src/Aligner.h	
${CMAKE_SOURCE_DIR}/src/KmerHistogram.h
${CMAKE_SOURCE_DIR}/src/Statistician.h
${CMAKE_SOURCE_DIR}/src/BestFirst.h
${CMAKE_SOURCE_DIR}/src/LockFreeQueue.h
${CMAKE_SOURCE_DIR}/src/AlignerParallel.h
${CMAKE_SOURCE_DIR}/src/IdentityCalculator.h
${CMAKE_SOURCE_DIR}/src/IdentityCalculator1.h
)

set(CLUSTER_HEADER_FILES 
${CMAKE_SOURCE_DIR}/src/meshclust/ClusteringUtil.h
${CMAKE_SOURCE_DIR}/src/meshclust/MeanShift.h
${CMAKE_SOURCE_DIR}/src/meshclust/MeanShiftLarge.h
${CMAKE_SOURCE_DIR}/src/meshclust/Cluster.h
${CMAKE_SOURCE_DIR}/src/meshclust/Reservoir.h
)



# RED source code
add_library(exception SHARED	
${CMAKE_SOURCE_DIR}/src/exception/FileDoesNotExistException.cpp  
${CMAKE_SOURCE_DIR}/src/exception/InvalidOperationException.cpp          
${CMAKE_SOURCE_DIR}/src/exception/InvalidScoreException.cpp
${CMAKE_SOURCE_DIR}/src/exception/InvalidInputException.cpp     
${CMAKE_SOURCE_DIR}/src/exception/InvalidOrderOfOperationsException.cpp  
${CMAKE_SOURCE_DIR}/src/exception/InvalidStateException.cpp
)

add_library(utility SHARED
${CMAKE_SOURCE_DIR}/src/utility/EmptyLocation.cpp	
${CMAKE_SOURCE_DIR}/src/utility/Location.cpp
${CMAKE_SOURCE_DIR}/src/utility/Util.cpp	
)
target_link_libraries(utility exception)

set(RED_HEADER_FILES  
${CMAKE_SOURCE_DIR}/src/utility/ILocation.h
${CMAKE_SOURCE_DIR}/src/nonltr/EnrichmentMarkovView.h
${CMAKE_SOURCE_DIR}/src/nonltr/KmerHashTable.h
${CMAKE_SOURCE_DIR}/src/nonltr/IChromosome.h
)

add_library(nonltr SHARED
${CMAKE_SOURCE_DIR}/src/nonltr/ChromosomeOneDigit.cpp
${CMAKE_SOURCE_DIR}/src/nonltr/HMM.cpp
${CMAKE_SOURCE_DIR}/src/nonltr/Scanner.cpp
${CMAKE_SOURCE_DIR}/src/nonltr/ChromDetectorMaxima.cpp
${CMAKE_SOURCE_DIR}/src/nonltr/ChromosomeRandom.cpp
${CMAKE_SOURCE_DIR}/src/nonltr/Scorer.cpp
${CMAKE_SOURCE_DIR}/src/nonltr/ChromListMaker.cpp
${CMAKE_SOURCE_DIR}/src/nonltr/DetectorMaxima.cpp
${CMAKE_SOURCE_DIR}/src/nonltr/LocationListCollection.cpp
${CMAKE_SOURCE_DIR}/src/nonltr/TableBuilder.cpp
${CMAKE_SOURCE_DIR}/src/nonltr/Chromosome.cpp
${CMAKE_SOURCE_DIR}/src/nonltr/LocationList.cpp
${CMAKE_SOURCE_DIR}/src/nonltr/Trainer.cpp
)

add_library(repeats SHARED
${CMAKE_SOURCE_DIR}/src/red/Red.cpp
)

add_executable(red ${CMAKE_SOURCE_DIR}/src/RepeatsDetector.cpp ${RED_HEADER_FILES})
target_link_libraries(red exception utility nonltr)

add_executable(redtest ${CMAKE_SOURCE_DIR}/src/red/RedTestDrive.cpp ${RED_HEADER_FILES})
target_link_libraries(redtest repeats exception utility nonltr)

add_executable(chromtest ${CMAKE_SOURCE_DIR}/src/test/chromosometest.cpp ${RED_HEADER_FILES})
target_link_libraries(chromtest nonltr utility)

# End of RED

add_executable(identity ${CMAKE_SOURCE_DIR}/src/Identity.cpp ${HEADER_FILES})
target_link_libraries(identity main)

add_executable(meshclust ${CMAKE_SOURCE_DIR}/src/meshclust/MeShClust.cpp ${HEADER_FILES} ${CLUSTER_HEADER_FILES})
target_link_libraries(meshclust clustering)

add_executable(makeclusters ${CMAKE_SOURCE_DIR}/src/apps/MakeClusters.cpp ${HEADER_FILES} ${CLUSTER_HEADER_FILES})
target_link_libraries(makeclusters main)

add_executable(generateclusters ${CMAKE_SOURCE_DIR}/src/apps/GenerateClusters.cpp ${HEADER_FILES})
target_link_libraries(generateclusters main)

##add_executable(identityclusteval ${CMAKE_SOURCE_DIR}/src/apps/IdentityClustEval.cpp ${HEADER_FILES})
##target_link_libraries(identityclusteval main)

add_executable(clusteval ${CMAKE_SOURCE_DIR}/src/apps/ClustEval.cpp ${HEADER_FILES} ${CLUSTER_HEADER_FILES})
target_link_libraries(clusteval clustering)

add_executable(clustevalalign ${CMAKE_SOURCE_DIR}/src/apps/ClustEvalAlignment.cpp ${HEADER_FILES} ${CLUSTER_HEADER_FILES})
target_link_libraries(clustevalalign clustering)

add_executable(mslToMatrix ${CMAKE_SOURCE_DIR}/src/msl/MslToMatrix.cpp ${MSL_SOURCE} ${HEADER_FILES})
target_link_libraries(mslToMatrix main)

add_executable(driver ${CMAKE_SOURCE_DIR}/src/test/Test.cpp ${TEST_SOURCE} ${HEADER_FILES})
target_link_libraries(driver main)

add_executable(ltrdetector ${CMAKE_SOURCE_DIR}/src/ltr/LTR.cpp ${HEADER_FILES} ${LTR_HEADER_FILES})
target_link_libraries(ltrdetector main ltr)

add_executable(look4ltrs ${CMAKE_SOURCE_DIR}/src/ltr/look4ltrs.cpp ${HEADER_FILES} ${LTR_HEADER_FILES})
target_link_libraries(look4ltrs main ltr repeats exception utility nonltr)

add_executable(look4ltrstest ${CMAKE_SOURCE_DIR}/src/ltr/MainTestDrive.cpp ${HEADER_FILES} ${LTR_HEADER_FILES})
target_link_libraries(look4ltrstest main ltr repeats exception utility nonltr)

#add_executable(identitytestcase ${CMAKE_SOURCE_DIR}/src/test/IdentityTestCase.cpp ${HEADER_FILES})
#target_link_libraries(identitytestcase main)