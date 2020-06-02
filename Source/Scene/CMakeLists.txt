#-----------------------------------------------------------------------------
# Create target
#-----------------------------------------------------------------------------
include(imstkAddLibrary)
imstk_add_library( Scene
  H_FILES
    imstkScene.h
    imstkCollisionGraph.h
	imstkcpdCollisionPair.h
	imstkcpdConstraintSet.h
	imstkcpdConstraintSolver.h
	imstkcpdExternalForceSolver.h
	imstkcpdScene.h
  CPP_FILES
    imstkScene.cpp
    imstkCollisionGraph.cpp
	imstkcpdCollisionPair.cpp
	imstkcpdConstraintSet.cpp
	imstkcpdConstraintSolver.cpp
	imstkcpdExternalForceSolver.cpp
	imstkcpdScene.cpp
  DEPENDS
    Core
    Collision
    SceneElements
    DynamicalModels
	Constraints
  )

#-----------------------------------------------------------------------------
# Testing
#-----------------------------------------------------------------------------
#if( ${PROJECT_NAME}_BUILD_TESTING )
#  add_subdirectory( Testing )
#endif()
