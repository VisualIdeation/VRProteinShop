/***********************************************************************
 VRProteinShop - Virtual reality version of ProtoShop (will hopefully be
 fully functional at some point).
 Copyright (c) 2002-2009 Oliver Kreylos

 This program is free software; you can redistribute it and/or modify it
 under the terms of the GNU General Public License as published by the
 Free Software Foundation; either version 2 of the License, or (at your
 option) any later version.

 This program is distributed in the hope that it will be useful, but
 WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 General Public License for more details.

 You should have received a copy of the GNU General Public License along
 with this program; if not, write to the Free Software Foundation, Inc.,
 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 ***********************************************************************/

#include <cstdio>
#include <cstring>
#include <string.h>
#include <iostream>
#include <stdexcept>

#include <Misc/ThrowStdErr.h>
#include <Misc/Timer.h>
#include <Comm/MulticastPipe.h>
#include <Math/Math.h>
#include <Geometry/Vector.h>
#include <Geometry/OrthonormalTransformation.h>
#include <Geometry/OrthogonalTransformation.h>
#include <GL/gl.h>
#include <GL/GLColorMap.h>
#include <GL/GLContextData.h>
#include <GL/GLLineIlluminator.h>
#include <GLMotif/Popup.h>
#include <GLMotif/Menu.h>
#include <GLMotif/SubMenu.h>
#include <GLMotif/PopupMenu.h>
#include <GLMotif/RowColumn.h>
#include <GLMotif/CascadeButton.h>
#include <GLMotif/TitleBar.h>
#include <GLMotif/PopupWindow.h>
#include <GLMotif/WidgetManager.h>
#include <Vrui/Vrui.h>
#include <Vrui/ClusterSupport.h>

#include <Protein/ClientServerPipe.h>
#include <Protein/CreateProtein.h>
#include <Protein/DragBox.h>
#include <Protein/EnergyAPI.h>
#include <Protein/Globals.h>
#include <Protein/ParsePdbFile.h>
#include <Protein/Protein.h>
#include <Protein/ProteinInteractor.h>
#include <Protein/ProteinRenderer.h>
#include <Protein/Stride.h>
#include <Protein/UndoBuffer.h>

#include "VRProteinShop.h"

/* Externs */

extern int offlineBuildBeta, checkBondSites, slowTrans;
extern double lastIKResidual;
extern float strand_to_coil_penalty, flatfrac, flatPhi, flflatPhi, raflatPhi,
   flatPsi, flflatPsi, raflatPsi, bb0[3], bb1[3];

/* Globals */

int alignBase = 2 ,planepr = 0;
const size_t inputfilename_size = 200;
char inputfilename[inputfilename_size];

/* Statics */

static int prolpr = 1;

/**********************************************
 Methods of class VRProteinShop::StructureLocator:
 **********************************************/

void VRProteinShop::StructureLocator::buttonPressCallback(
		Vrui::LocatorTool::ButtonPressCallbackData* cbData) {
	MD::Point p = cbData->currentTransformation.getOrigin();

	{
		Threads::Mutex::Lock proteinLock(application->proteinMutex);
		switch (application->selectionMode) {
		case SELECT_SECONDARYSTRUCTURE: {
			/* Find the secondary structure containing the current locator position: */
			MD::Protein::StructureSelector ss =
					application->protein->pickStructure(p);

			/* Select the picked secondary structure: */
			application->interactor->selectStructure(ss);
			break;
		}

		case SELECT_SINGLERESIDUE: {
			/* Find the residue containing the current locator position: */
			//TODO:			MD::Protein::Residue* res=application->protein->pickResidue(p);

			/* Select the picked residue: */
			//TODO:			application->interactor->setSelectedResidue(res);
			break;
		}

		case SELECT_TOGGLEACTIVECOIL: {
			/* Find the secondary structure containing the current locator position: */
			MD::Protein::StructureSelector ss =
					application->protein->pickStructure(p);

			/* Toggle the activation status of the picked secondary structure: */
			application->interactor->toggleCoil(ss);
			break;
		}
		}

		/* Invalidate the graphics cache: */
		++application->graphicsVersion;
	}
}

/**********************************************
 Methods of class VRProteinShop::StructureDragger:
 **********************************************/

void VRProteinShop::StructureDragger::dragStartCallback(
		Vrui::DraggingTool::DragStartCallbackData* cbData) {
	/* Try picking the structure dragging box: */
	isDragging = false;
	if (cbData->rayBased) {
		rayBased = true;
		initialRay = cbData->ray;
		initialRay.inverseTransform(cbData->startTransformation);

		DragBox::Point start(cbData->ray.getOrigin());
		DragBox::Vector direction(cbData->ray.getDirection());
		//TODO:		isDragging=application->interactor->pickDragBox(start,direction);
	} else {
		rayBased = false;

		DragBox::OTransformation dt(
				cbData->startTransformation.getTranslation(),
				cbData->startTransformation.getRotation());
		isDragging = application->interactor->pickDragBox(dt);
	}
	if (isDragging) {
		/* Prepare the IK thread: */
		{
			Threads::Mutex::Lock proteinLock(application->proteinMutex);
			application->ikUpdateRequested = false;
			application->ikUpdatePosted = false;
		}
	}
}

void VRProteinShop::StructureDragger::dragCallback(
		Vrui::DraggingTool::DragCallbackData* cbData) {
	if (isDragging) {
		/* Drag the structure drag box: */
		{
			Threads::Mutex::Lock proteinLock(application->proteinMutex);
			if (rayBased) {
				Vrui::Ray ray = initialRay;
				ray.transform(cbData->currentTransformation);
				DragBox::Point start(ray.getOrigin());
				DragBox::Vector direction(ray.getDirection());
				//TODO:			application->interactor->dragBox(start,direction);
			} else {
				DragBox::OTransformation dt(
						cbData->currentTransformation.getTranslation(),
						cbData->currentTransformation.getRotation());
				application->interactor->dragBox(dt);
			}

			/* Post an update request to the IK thread: */
			application->ikUpdateRequested = true;
			if (application->ikUpdateDone)
				application->ikUpdateRequestedCond.signal();
		}
	}
}

void VRProteinShop::StructureDragger::dragEndCallback(
		Vrui::DraggingTool::DragEndCallbackData* cbData) {
	if (isDragging) {
		isDragging = false;

		/* Stop the IK thread: */
		{
			Threads::Mutex::Lock proteinLock(application->proteinMutex);
			application->ikUpdateRequested = false;
		}
	}
}

/****************************
 Methods of class VRProteinShop:
 ****************************/

void VRProteinShop::addProtein (MD::Protein* newProtein, int Id)
{
    // called to add a new protein to the workspace
    using namespace MD;

    if ( !newProtein ) return;
    ProteinState * state = createProtein(inputfilename);
    if ( !state ) return;

    state->protein = newProtein;
    configFile->setCurrentSection ("/ProteinRenderer");
    state->proteinRenderer = new ProteinRenderer(configFile->getCurrentSection(), state->protein);
    configFile->setCurrentSection("/");

    state->proteinId = Id;
	if ( validPdbForEnergy )
    {
        msg.Debug (PDBID, "This PDB file has the correct number of atoms.");
        if ( energyLibrary )
        {
            state->energyCalculator = energyLibrary->createEnergyCalculator (state->protein);
            if ( state->energyCalculator )
            {
                state->visualizeEnergy = false;
                state->visualizeEnergyMinRange = 0.0;
                state->visualizeEnergyMaxRange = 1.0;
                state->energyRenderer = new EnergyRenderer (state);
            }
            else
            {
                msg.Warn (PDBID, "There is no Energy Calculator.");
            }
        }
    }
    else
    {
        // some hydrogen atoms may be missing from this molecule

    }
    state->proteinRenderer->setMapAtomValues (state->visualizeEnergy);
    if ( state->visualizeEnergy )
    {
        state->proteinRenderer->setMapAtomValueRange (
            state->visualizeEnergyMinRange,
            state->visualizeEnergyMaxRange
        );
    }
    state->interactor = new ProteinInteractor (
        state->protein,
        state->proteinRenderer,
        state->energyCalculator,
        undoBuffer
    );
}

GLMotif::Popup* VRProteinShop::createGlobalDisplayMenu(void) {
	GLMotif::Popup* globalDisplayMenuPopup = new GLMotif::Popup(
			"GlobalDisplayMenuPopup", Vrui::getWidgetManager());

	GLMotif::SubMenu* globalDisplayMenu = new GLMotif::SubMenu("DisplayMenu",
			globalDisplayMenuPopup, false);

	globalDrawAtomsToggle = new GLMotif::ToggleButton("AtomsToggle",
			globalDisplayMenu, "Draw atoms");
	globalDrawAtomsToggle->getValueChangedCallbacks().add(this,
			&VRProteinShop::menuToggleSelectCallback);
	globalDrawBondsToggle = new GLMotif::ToggleButton("BondsToggle",
			globalDisplayMenu, "Draw sidechain bonds");
	globalDrawBondsToggle->getValueChangedCallbacks().add(this,
			&VRProteinShop::menuToggleSelectCallback);
	globalDrawBackboneToggle = new GLMotif::ToggleButton("BackboneToggle",
			globalDisplayMenu, "Draw backbone ribbon");
	globalDrawBackboneToggle->getValueChangedCallbacks().add(this,
			&VRProteinShop::menuToggleSelectCallback);
	globalDrawCartoonToggle = new GLMotif::ToggleButton("CartoonToggle",
			globalDisplayMenu, "Draw structure cartoons");
	globalDrawCartoonToggle->getValueChangedCallbacks().add(this,
			&VRProteinShop::menuToggleSelectCallback);
	globalDrawHydrogenBondsToggle = new GLMotif::ToggleButton(
			"HydrogenBondsToggle", globalDisplayMenu, "Draw hydrogen bonds");
	globalDrawHydrogenBondsToggle->getValueChangedCallbacks().add(this,
			&VRProteinShop::menuToggleSelectCallback);
	globalDrawHydrogenBondSitesToggle = new GLMotif::ToggleButton(
			"HydrogenBondSitesToggle", globalDisplayMenu,
			"Draw hydrogen bond sites");
	globalDrawHydrogenBondSitesToggle->getValueChangedCallbacks().add(this,
			&VRProteinShop::menuToggleSelectCallback);
	globalDrawHydrogenCagesToggle = new GLMotif::ToggleButton(
			"HydrogenCagesToggle", globalDisplayMenu, "Draw hydrogen cages");
	globalDrawHydrogenCagesToggle->getValueChangedCallbacks().add(this,
			&VRProteinShop::menuToggleSelectCallback);
	globalDrawAtomCollisionsToggle = new GLMotif::ToggleButton(
			"AtomCollisionsToggle", globalDisplayMenu,
			"Visualize atom collisions");
	globalDrawAtomCollisionsToggle->getValueChangedCallbacks().add(this,
			&VRProteinShop::menuToggleSelectCallback);
	if (densityRenderer != 0) {
		globalDrawDensityToggle = new GLMotif::ToggleButton("DensityToggle",
				globalDisplayMenu, "Visualize density distribution");
		globalDrawDensityToggle->getValueChangedCallbacks().add(this,
				&VRProteinShop::menuToggleSelectCallback);
	} else
		globalDrawDensityToggle = 0;

	globalDisplayMenu->manageChild();

	return globalDisplayMenuPopup;
}

GLMotif::PopupMenu* VRProteinShop::createMainMenu(void) {
	GLMotif::PopupMenu* mainMenuPopup = new GLMotif::PopupMenu("MainMenuPopup",
			Vrui::getWidgetManager());
	mainMenuPopup->setTitle("GLMotif:: ProtoShop");

	GLMotif::Menu* mainMenu = new GLMotif::Menu("MainMenu", mainMenuPopup,
			false);

	GLMotif::Button* centerDisplayButton = new GLMotif::Button(
			"CenterDisplayButton", mainMenu, "Center Display");
	centerDisplayButton->getSelectCallbacks().add(this,
			&VRProteinShop::centerDisplayCallback);

	GLMotif::CascadeButton* globalDisplayMenuCascade =
			new GLMotif::CascadeButton("GlobalDisplayMenuCascade", mainMenu,
					"Global Display Settings");
	globalDisplayMenuCascade->setPopup(createGlobalDisplayMenu());

	GLMotif::ToggleButton* showInteractionDialogToggle =
			new GLMotif::ToggleButton("ShowInteractionDialogToggle", mainMenu,
					"Show Interaction Dialog");
	showInteractionDialogToggle->setToggle(false);
	showInteractionDialogToggle->getValueChangedCallbacks().add(this,
			&VRProteinShop::showInteractionDialogValueChangedCallback);

	mainMenu->manageChild();

	return mainMenuPopup;
}

GLMotif::PopupWindow* VRProteinShop::createInteractionDialog(void) {
	GLMotif::PopupWindow* interactionPopup = new GLMotif::PopupWindow(
			"InteractionPopup", Vrui::getWidgetManager(), "Interaction Dialog");

	GLMotif::RowColumn* interactionPanel = new GLMotif::RowColumn(
			"InteractionPanel", interactionPopup, false);

	GLMotif::RadioBox* selectionModeBox = new GLMotif::RadioBox(
			"SelectionModeBox", interactionPanel, false);
	selectionModeBox->setSelectionMode(GLMotif::RadioBox::ALWAYS_ONE);

	selectionModeBox->addToggle("Secondary Structure");
	selectionModeBox->addToggle("Single Residue");
	selectionModeBox->addToggle("Toggle Active Coil Region");

	selectionModeBox->manageChild();
	switch (selectionMode) {
	case SELECT_SECONDARYSTRUCTURE:
		selectionModeBox->setSelectedToggle(0);
		break;

	case SELECT_SINGLERESIDUE:
		selectionModeBox->setSelectedToggle(1);
		break;

	case SELECT_TOGGLEACTIVECOIL:
		selectionModeBox->setSelectedToggle(2);
		break;
	}
	selectionModeBox->getValueChangedCallbacks().add(this,
			&VRProteinShop::selectionModeValueChangedCallback);

	GLMotif::RadioBox* ikUpdateDirectionBox = new GLMotif::RadioBox(
			"IkUpdateDirectionBox", interactionPanel, false);
	ikUpdateDirectionBox->setSelectionMode(GLMotif::RadioBox::ALWAYS_ONE);

	ikUpdateDirectionBox->addToggle("Left (C to N)");
	ikUpdateDirectionBox->addToggle("Right (N to C)");

	ikUpdateDirectionBox->manageChild();
	switch (updateDirection) {
	case UPDATE_LEFT:
		ikUpdateDirectionBox->setSelectedToggle(0);
		break;

	case UPDATE_RIGHT:
		ikUpdateDirectionBox->setSelectedToggle(0);
		break;
	}
	ikUpdateDirectionBox->getValueChangedCallbacks().add(this,
			&VRProteinShop::ikUpdateDirectionValueChangedCallback);

	GLMotif::Button* resetActiveCoilRegionsButton = new GLMotif::Button(
			"ResetActiveCoilRegionsButton", interactionPanel,
			"Reset Active Coil Regions");
	resetActiveCoilRegionsButton->getSelectCallbacks().add(this,
			&VRProteinShop::resetActiveCoilRegionsCallback);

	autoBondingBox = new GLMotif::RadioBox("AutoBondingBox", interactionPanel,
			false);
	autoBondingBox->setSelectionMode(GLMotif::RadioBox::ATMOST_ONE);

	autoBondingBox->addToggle("Parallel");
	autoBondingBox->addToggle("Anti-parallel");

	autoBondingBox->manageChild();
	autoBondingBox->setSelectedToggle(0);
	autoBondingBox->getValueChangedCallbacks().add(this,
			&VRProteinShop::autoBondingValueChangedCallback);

	interactionPanel->manageChild();

	return interactionPopup;
}

void VRProteinShop::updateMenuState(void) {
	globalDrawAtomsToggle->setToggle(proteinRenderer->getDrawAtoms());
	globalDrawBondsToggle->setToggle(proteinRenderer->getDrawBonds());
	globalDrawBackboneToggle->setToggle(
			proteinRenderer->getDrawBackboneRibbon());
	globalDrawCartoonToggle->setToggle(proteinRenderer->getDrawCartoon());
	globalDrawHydrogenBondsToggle->setToggle(
			proteinRenderer->getDrawHydrogenBonds());
	globalDrawHydrogenBondSitesToggle->setToggle(
			proteinRenderer->getDrawHydrogenBondSites());
	globalDrawHydrogenCagesToggle->setToggle(
			proteinRenderer->getDrawHydrogenCages());
	globalDrawAtomCollisionsToggle->setToggle(
			proteinRenderer->getDrawCollisions());
	if (densityRenderer != 0)
		globalDrawDensityToggle->setToggle(drawDensity);
}

void* VRProteinShop::ikUpdateThreadMethod(void) {
	/* Enable immediate cancellation of this thread: */
	Threads::Thread::setCancelState(Threads::Thread::CANCEL_ENABLE);
	Threads::Thread::setCancelType(Threads::Thread::CANCEL_ASYNCHRONOUS);

	ikUpdateDone = true;
	proteinMutex.lock();
	while (true) {
		/* Wait for the next update request: */
		if (!ikUpdateRequested) {
			ikUpdateDone = true;
			ikUpdateDoneCond.signal();
			ikUpdateRequestedCond.wait(proteinMutex);
		}

		/* Grab the current goal transformation from the drag box: */
		ProteinInteractor::Transformation goalTransformation =
				interactor->getDragTransformation();
		// ikUpdateRequested=false; // Keep going even if the goal transformation did not change
		ikUpdateDone = false;
		proteinMutex.unlock();

		/* Perform IK steps: */
		interactor->drag(goalTransformation);

		/* Apply changes to protein: */
		proteinMutex.lock();
		interactor->applyChanges(); // Notifies server of changes if clientPipe!=0
		++graphicsVersion;
		ikUpdatePosted = true;
	}

	return 0;
}

VRProteinShop::VRProteinShop(int& argc, char**& argv, char**& appDefaults) :
	Vrui::Application(argc, argv, appDefaults), protein(0), proteinID(0),
			graphicsVersion(1), proteinRenderer(0), drawDensity(false),
			densityRenderer(0), densityPalette(0),
			selectionMode(SELECT_SECONDARYSTRUCTURE), autoBondingMode(AUTOBOND_NONE),
			updateDirection(UPDATE_RIGHT), undoBuffer(0), interactor(0),
			interactorPipe(Vrui::openPipe()), mainMenu(0), interactionDialog(0), validPdbForEnergy(false) {
	/* Scan the command line: */
	const char* inputFileName = 0;
	const char* densityFileName = 0;
	const char* paletteFileName = 0;
	for (int i = 1; i < argc; ++i) {
		if (argv[i][0] == '-') {
			if (strcasecmp(argv[i] + 1, "vol") == 0 || strcasecmp(argv[i] + 1,
					"fvol") == 0) {
				/* Treat the next argument as a density file name: */
				++i;
				densityFileName = argv[i];
			} else if (strcasecmp(argv[i] + 1, "pal") == 0) {
				/* Treat the next argument as a palette file name: */
				++i;
				paletteFileName = argv[i];
			}
		} else {
			/* Treat this argument as a protein file name: */
			inputFileName = argv[i];
		}
	}
	if (inputFileName == 0)
		throw std::runtime_error(
				"No PDB file name or prediction file name provided");

    try
    {
        /* Open configuration file: */
        configFile = new Misc::ConfigurationFile(VRPROTEINSHOP_CONFIGFILENAME);
    }
    catch ( std::runtime_error error )
    {
        msg.Error (CONFIGID,"",error.what());
        std::cerr << "Could not load configuration file ProteinShop.cfg" << "\n";
    }

	/* Query extension of input file: */
	const char* extPtr;
	for (const char* cPtr = inputFileName; *cPtr != '\0'; ++cPtr)
		if (*cPtr == '.')
			extPtr = cPtr;
	if (strcasecmp(extPtr, ".pdb") == 0) {
		/* Load a protein model: */
		protein = MD::parsePdbFile(inputFileName);
	} else if (strcasecmp(extPtr, ".pred") == 0) {
		/* Build a protein model: */
		MD::ReadStandards( VRPROTEINSHOP_STANDARDSDIR);
		protein = MD::ReadPredictionFile(inputFileName, 1);
	} else
		throw std::runtime_error(
				std::string("Unrecognized file name extension") + std::string(
						extPtr));

	/* Initialize dependent protein state: */
	proteinRenderer = new MD::ProteinRenderer(configFile->getSection("/ProteinRenderer"), protein);
	energyCalculator = 0;
	undoBuffer = new UndoBuffer;
	interactor = new ProteinInteractor(protein, proteinRenderer, energyCalculator, undoBuffer);

#if 0
	/* Disable bond cylinder rendering for all structures, but enable it for the entire protein: */
	for(int i=0;i<protein->getNumStructures();++i)
	proteinRenderer->setDrawBonds(protein->pickStructure(i),false);
	proteinRenderer->setDrawBonds(true);
#endif

	/* Initialize inverse kinematics state: */
	ikUpdateRequested = false;
	ikUpdateDone = false;
	ikUpdatePosted = false;
	ikUpdateThread.start(this, &VRProteinShop::ikUpdateThreadMethod);

	if (densityFileName != 0) {
		/* Load palette renderer from file: */
		densityRenderer = new PaletteRenderer();

		/* Initialize renderer settings: */
		// densityRenderer->setUseNPOTDTextures(true);
		densityRenderer->setVoxelAlignment(VolumeRenderer::VERTEX_CENTERED);
		// densityRenderer->setRenderingMode(VolumeRenderer::AXIS_ALIGNED);
		densityRenderer->setRenderingMode(VolumeRenderer::VIEW_PERPENDICULAR);
		densityRenderer->setInterpolationMode(VolumeRenderer::LINEAR);
		densityRenderer->setTextureFunction(VolumeRenderer::REPLACE);
		densityRenderer->setSliceFactor(VolumeRenderer::Scalar(1.4142));
		densityRenderer->setAutosaveGLState(true);
		densityRenderer->setTextureCaching(true);

		/* Create a color map: */
		if (paletteFileName != 0)
			densityPalette = new GLColorMap(paletteFileName, 0.0f, 1.0f);
		else
			densityPalette = new GLColorMap(GLColorMap::RAINBOW
					| GLColorMap::RAMP_ALPHA, 1.0f, 1.0f, 0.0f, 1.0f);
		densityPalette->changeTransparency(1.0);
		densityPalette->premultiplyAlpha();
		densityRenderer->setSharePalette(false);
		densityRenderer->setColorMap(densityPalette);
	}

	/* Create the main menu: */
	mainMenu = createMainMenu();
	Vrui::setMainMenu(mainMenu);
	updateMenuState();

	/* Create the interaction dialog: */
	interactionDialog = createInteractionDialog();
	// Vrui::popupPrimaryScreenWidget(interactionDialog,0.1,0.1);
	// Vrui::getWidgetManager()->hide(interactionDialog);

	/* Initialize the navigation transformation: */
	centerDisplayCallback(0);
}

VRProteinShop::~VRProteinShop(void) {
	delete mainMenu;
	delete interactionDialog;

	/* Kill the IK update thread: */
	ikUpdateThread.cancel();
	ikUpdateThread.join();

	delete interactor;
	delete interactorPipe;
	delete undoBuffer;
	delete densityPalette;
	delete densityRenderer;
	delete proteinRenderer;
	delete protein;
}

void VRProteinShop::initContext(GLContextData& contextData) const {
	/* Create a new context entry: */
	contextData.addDataItem(this, new DataItem);
}

void VRProteinShop::loadProteins (const char *filename)
{
    using namespace MD;
    Protein *protein = 0;
    Stride* pred = new Stride;
	char* chain = new char[1];
	int models = parseModels(filename);
	const char* chainList = parseChains(filename);
	char* selectedModel = new char[models];
	char* selectedChain = new char[strlen(chainList)];

    if ( filename )
    {
		if(!pred->stridePrediction(filename))
		{
			std::cout << "Not able to predict the second structure!" << "\n";
			return;
		}

		strncpy (inputfilename, filename, inputfilename_size);
        inputfilename[inputfilename_size-1] = 0;
		try
        {
            /* Check the requirements of energy computation for every PDB file: */
            validPdbForEnergy = checkPdbFile (filename);
        }
        catch ( std::runtime_error error )
        {
            /* Silently disable the EnergyDialog */
            validPdbForEnergy = false;
            std::cout << "PDB file not valid for energy computation." << "\n";
        }

		if ( energyLibrary && !validPdbForEnergy)
            msg.Warn (PDBID, "This PDB file is not ready for Energy Calculation.");

		try
        {
   			// chain selection
			if(strlen(chainList)>1)
			{
			msg.Info(PDBID, "This is a X-ray diffraction solved structure. ", filename);
			strcpy(selectedChain,configFile->retrieveString("/ProteinCreator/chainSelection","").c_str());
			if (!strlen(selectedChain))
			{
				// all the chains will be loaded
				for(int j =0; j< strlen(chainList); j++)
				{
					strncpy(chain, &chainList[j], 1);
					chain[1]='\0';
					// Load a protein model: Disable ACE/NME cap if this pdb file is not valid
           			if(validPdbForEnergy)
						protein = parsePdbFile (filename, chain);
					else
						protein = parsePdbFile ( filename, false, chain);
					if ( protein ) addProtein (protein, proteinID++);
				}
			}
			else
			{
				selectedChain = strtok( selectedChain, " " );
				while( selectedChain != NULL )
				{
					if (sscanf( selectedChain, "%c", chain ))
						chain[1]='\0';

					selectedChain = strtok( NULL, " " );
			    	bool found = false;
					// load selected chain
					for(int j =0; j< strlen(chainList); j++)
        			{
						// make sure the existness of selected chain
						if(strncmp(chain, &chainList[j], 1)==0)
						{
            				if(validPdbForEnergy)
								protein = parsePdbFile (filename, chain);
							else
								protein = parsePdbFile ( filename, false, chain);
							if ( protein ) addProtein (protein, proteinID++);
							found = true;
						}
					}
					if(!found)
						msg.Warn(PDBID, "This chain (", chain, ") was not found! Check configuration file.");
				}
			}
			}
			else // model selection
			{
				msg.Info(PDBID, "The structure was determined using NMR spectroscopy. ", filename);
				strcpy(selectedModel,configFile->retrieveString("/ProteinCreator/multipleModel","").c_str());
				if ( !strlen(selectedModel) )
           		{
					// all the models will be loaded
					for (int i = 1; i< models+1; i++)
					{
						if(validPdbForEnergy)
							protein = parsePdbFile (filename, i);
						else
							protein = parsePdbFile ( filename, false, 0, i);
						if ( protein ) addProtein (protein, proteinID++);
 					}
				}
				else
				{
 					int number;
					selectedModel = strtok( selectedModel, " " );
					while( selectedModel != NULL )
					{
						if (!sscanf(selectedModel , "%d", &number ))
							msg.Warn(PDBID, "No number found!");

            	  		selectedModel = strtok( NULL, " " );

						if(number <0 || number > models)
						{
							sprintf(msg.buf, "This model (%d) was not found! Check configuration file.", number);
							msg.Warn(PDBID, msg.buf);
							continue;
						}
						if(validPdbForEnergy)
							protein = parsePdbFile (filename, number);
						else
							protein = parsePdbFile ( filename, false, 0, number);
						if ( protein ) addProtein (protein, proteinID++);
					}
				}
			}
		}
        catch ( std::runtime_error error )
        {
            std::cerr << "Error loading PDB file:\n"<< error.what() << "\n";
            msg.Error (PDBID, error.what());
       }
    }
	delete pred;
	delete [] chain;
	delete [] selectedModel;
	delete [] selectedChain;
}

MD::Protein * VRProteinShop::loadPdb (const char *filename)
{
    using namespace MD;
    Protein *protein = 0;
    Stride* pred = new Stride;
    if ( filename )
    {
		if(!pred->stridePrediction(filename))
		{
			std::cout << "Not able to predict the second structure!" << "\n";
			return 0;
		}
        strncpy (inputfilename, filename, inputfilename_size);
        inputfilename[inputfilename_size-1] = 0;
        try
        {
            /* Check the requirements of energy computation for every PDB file: */
            validPdbForEnergy = checkPdbFile (filename);
        }
        catch ( std::runtime_error error )
        {
            /* Silently disable the EnergyDialog */
            validPdbForEnergy = false;
            std::cout << "PDB file not valid for energy computation." << "\n";
        }
		if ( energyLibrary && !validPdbForEnergy)
            msg.Warn (PDBID, "This PDB file is not ready for Energy Calculation.");
        try
        {
            /* Load a protein model: */
			if(validPdbForEnergy)
				protein = parsePdbFile (filename);
			else
				protein = parsePdbFile (filename, false);
        }
        catch ( std::runtime_error error )
        {
            std::cerr << "Error loading PDB file:\n%" << error.what() << "\n";
            msg.Error (PDBID, error.what());
        }
    }
	delete pred;
    return protein;
}

void VRProteinShop::toolCreationCallback(
		Vrui::ToolManager::ToolCreationCallbackData* cbData) {
	/* Check if the new tool is a locator tool: */
	Vrui::LocatorTool* locatorTool =
			dynamic_cast<Vrui::LocatorTool*> (cbData->tool);
	if (locatorTool != 0) {
		/* Create a structure locator object and associate it with the new tool: */
		StructureLocator* newLocator = new StructureLocator(locatorTool, this);

		/* Add new locator to list: */
		locators.push_back(newLocator);
	}

	/* Check if the new tool is a dragging tool: */
	Vrui::DraggingTool* draggingTool =
			dynamic_cast<Vrui::DraggingTool*> (cbData->tool);
	if (draggingTool != 0) {
		/* Create a structure dragger object and associate it with the new tool: */
		StructureDragger* newDragger = new StructureDragger(draggingTool, this);

		/* Add new dragger to list: */
		draggers.push_back(newDragger);
	}
}

void VRProteinShop::toolDestructionCallback(
		Vrui::ToolManager::ToolDestructionCallbackData* cbData) {
	/* Check if the to-be-destroyed tool is a locator tool: */
	Vrui::LocatorTool* locatorTool =
			dynamic_cast<Vrui::LocatorTool*> (cbData->tool);
	if (locatorTool != 0) {
		/* Find the locator associated with the tool in the list: */
		LocatorList::iterator lIt;
		for (lIt = locators.begin(); lIt != locators.end(); ++lIt)
			if ((*lIt)->getTool() == locatorTool) {
				/* Remove the locator: */
				delete *lIt;
				locators.erase(lIt);
				break;
			}
	}

	/* Check if the to-be-destroyed tool is a dragging tool: */
	Vrui::DraggingTool* draggingTool =
			dynamic_cast<Vrui::DraggingTool*> (cbData->tool);
	if (draggingTool != 0) {
		/* Find the dragger associated with the tool in the list: */
		DraggerList::iterator dIt;
		for (dIt = draggers.begin(); dIt != draggers.end(); ++dIt)
			if ((*dIt)->getTool() == draggingTool) {
				/* Remove the dragger: */
				delete *dIt;
				draggers.erase(dIt);
				break;
			}
	}
}

void VRProteinShop::frame(void) {
	if (ikUpdatePosted) {
		ikUpdatePosted = false;

		if (ikUpdateDone) {
			/* Snap the drag box back to its true position: */
			interactor->releaseDragBox();
			interactor->finishInteraction();
		}
	}

	/* Set the viewing direction for illuminated lines and density rendering: */
	viewDirection = VolumeRenderer::Vector(Vrui::getViewDirection());
	GLLineIlluminator::Vector lineViewDirection(viewDirection.getComponents());
	proteinRenderer->setViewDirection(lineViewDirection);
	proteinRenderer->setLightDirection(lineViewDirection);
}

void VRProteinShop::display(GLContextData& contextData) const {
	/* Retrieve context entry: */
	DataItem* dataItem = contextData.retrieveDataItem<DataItem> (this);

	/* Draw the protein: */
	{
		Threads::Mutex::Lock proteinLock(proteinMutex);
#if 1
		proteinRenderer->glRenderAction(contextData);
#else
		if(dataItem->graphicsVersion!=graphicsVersion)
		{
			glNewList(dataItem->graphicsDisplayListId,GL_COMPILE_AND_EXECUTE);
			proteinRenderer->glRenderAction(contextData);
			glEndList();
			dataItem->graphicsVersion=graphicsVersion;
		}
		else
		glCallList(dataItem->graphicsDisplayListId);
#endif
	}

	if (densityRenderer != 0 && drawDensity) {
		/* Render the density distribution: */
		densityRenderer->renderBlock(contextData, viewDirection);
	}

	/* Draw the structure selection box: */
	glPushAttrib(GL_LIGHTING_BIT);
	glDisable(GL_LIGHTING);
	interactor->glRenderAction(contextData);
	glPopAttrib();
}

void VRProteinShop::centerDisplayCallback(Misc::CallbackData*) {
	/* Calculate the model's position and size: */
	Vrui::Point modelCenter = Vrui::Point(protein->calcCentroid());
	Vrui::Scalar modelSize = Vrui::Scalar(protein->calcRadius());

	/* Set the navigation transformation: */
	Vrui::setNavigationTransformation(modelCenter, modelSize);
}

void VRProteinShop::menuToggleSelectCallback(
		GLMotif::ToggleButton::ValueChangedCallbackData* cbData) {
	if (cbData->toggle == globalDrawAtomsToggle)
		proteinRenderer->setDrawAtoms(cbData->toggle->getToggle());
	else if (cbData->toggle == globalDrawBondsToggle)
		proteinRenderer->setDrawBonds(cbData->toggle->getToggle());
	else if (cbData->toggle == globalDrawBackboneToggle)
		proteinRenderer->setDrawBackboneRibbon(cbData->toggle->getToggle());
	else if (cbData->toggle == globalDrawCartoonToggle) {
		proteinRenderer->setDrawBackbone(!cbData->toggle->getToggle());
		proteinRenderer->setDrawCartoon(cbData->toggle->getToggle());
	} else if (cbData->toggle == globalDrawHydrogenBondsToggle)
		proteinRenderer->setDrawHydrogenBonds(cbData->toggle->getToggle());
	else if (cbData->toggle == globalDrawHydrogenBondSitesToggle)
		proteinRenderer->setDrawHydrogenBondSites(cbData->toggle->getToggle());
	else if (cbData->toggle == globalDrawHydrogenCagesToggle)
		proteinRenderer->setDrawHydrogenCages(cbData->toggle->getToggle());
	else if (cbData->toggle == globalDrawAtomCollisionsToggle)
		proteinRenderer->setDrawCollisions(cbData->toggle->getToggle());
	else if (cbData->toggle == globalDrawDensityToggle)
		drawDensity = cbData->toggle->getToggle();

	updateMenuState();

	{
		/* Invalidate the graphics cache: */
		Threads::Mutex::Lock proteinLock(proteinMutex);
		++graphicsVersion;
	}
}

void VRProteinShop::showInteractionDialogValueChangedCallback(
		GLMotif::ToggleButton::ValueChangedCallbackData* cbData) {
	if (cbData->set) {
		/* Pop up the interaction dialog at the same position as the main menu: */
		Vrui::getWidgetManager()->popupPrimaryWidget(interactionDialog,
				Vrui::getWidgetManager()->calcWidgetTransformation(mainMenu));
	} else {
		Vrui::popdownPrimaryWidget(interactionDialog);
	}
}

void VRProteinShop::selectionModeValueChangedCallback(
		GLMotif::RadioBox::ValueChangedCallbackData* cbData) {
	switch (cbData->radioBox->getToggleIndex(cbData->newSelectedToggle)) {
	case 0:
		selectionMode = SELECT_SECONDARYSTRUCTURE;
		break;

	case 1:
		selectionMode = SELECT_SINGLERESIDUE;
		break;

	case 2:
		selectionMode = SELECT_TOGGLEACTIVECOIL;
		break;
	}
}

void VRProteinShop::ikUpdateDirectionValueChangedCallback(
		GLMotif::RadioBox::ValueChangedCallbackData* cbData) {
}

void VRProteinShop::resetActiveCoilRegionsCallback(Misc::CallbackData*) {
}

void VRProteinShop::autoBondingValueChangedCallback(
		GLMotif::RadioBox::ValueChangedCallbackData* cbData) {
}

/****************
 Global functions:
 ****************/

void updateDihedralAngles(void)
{
}

void updateProteinNow(void) {
}

void applyTransformation(const ProteinInteractor::Transformation& transformation)
{
    int rotatePrint = 0;
    using namespace MD;
    ProteinState *state = curProtein();
    if ( !state ) {
       printf("applyTransormation returning early because of no state\n");
       return;
       }

    typedef ProteinInteractor::Transformation Transformation;
    typedef Transformation::Rotation Rotation;
    typedef Transformation::Point Point;
    typedef Transformation::Scalar Scalar;
    typedef Transformation::Vector Vector;

    /* Decompose the given transformation: */
    Vector translation=transformation.getTranslation();
    Vector rotationAxis=transformation.getRotation().getAxis();
    Scalar rotationAngle=transformation.getRotation().getAngle();
    Point rotateCenterStart=state->interactor->getBoxRotateCenter();
    Point rotateCenterEnd=transformation.transform(rotateCenterStart);

    /* Calculate partial transformations to go in several small steps: */
    // The values of 1 below were previously 2, but the change should help prevent flipping
    int numSteps=int(Math::ceil(Geometry::dist(rotateCenterStart,rotateCenterEnd)/1.0)); // At most 1 Angstrom per step
    // int numRotateSteps=int(Math::ceil(Math::abs(rotationAngle)/0.1)); // At most 0.1 radians per step
    int numRotateSteps=int(Math::ceil(Math::abs(rotationAngle)/0.05)); // At most 0.05 radians per step
    if(rotatePrint) printf("numSteps %d numRotateSteps %d\n",
       numSteps, numRotateSteps);
    if(numSteps<numRotateSteps)
        numSteps=numRotateSteps;
    if(state->interactor->startInteraction())
    {
        for(int i=0;i<numSteps;++i)
        {
            Scalar s=Scalar(i+1)/Scalar(numSteps);
            Rotation r=Rotation::rotateAxis(rotationAxis,rotationAngle*s);
            Point rc=Geometry::affineCombination(rotateCenterStart,rotateCenterEnd,s);
            Vector t=rc-r.transform(rotateCenterStart);
            Transformation stepTransformation(t,r);
            state->interactor->drag(stepTransformation);
            state->interactor->applyChanges();
	    if(0) if(checkBondAngles(7))
	       printf("bond angle changed after step %d of %d in applyTransformation\n",
	          i, numSteps);
            if(offlineBuildBeta == 0) updateProteinNow();
	    if(rotatePrint) {
	       printf("%d ", i);
	       if(slowTrans) {
	          printf("%g lastIKResidual\n", lastIKResidual);
		  if(i >= numSteps - 3) {
	             fflush(stdout);
		     getchar();
		     }
		  }
	       }
        }
	if(rotatePrint) printf("\n");
        state->interactor->finishInteraction();
        if(offlineBuildBeta == 0) updateDihedralAngles();
    }
    else
        printf("state->interactor->startInteraction() returned 0 in applyTransformation.\n");
}

int zipAntiParallel(MD::Protein::Residue* manipulatedResidue,MD::Protein::Residue* anchorResidue, float pl[3], float cn[3])
{
    typedef ProteinInteractor::Transformation Transformation;
    float test, dentest;

    /* Zip the two selected residues together: */
    MD::Protein::Dipole amide1=manipulatedResidue->getAmide();
    MD::Protein::Dipole carboxyl1=manipulatedResidue->getCarboxyl();
    MD::Protein::Dipole amide2=anchorResidue->getAmide();
    MD::Protein::Dipole carboxyl2=anchorResidue->getCarboxyl();
    if(amide1.isValid()&&carboxyl1.isValid()&&amide2.isValid()&&carboxyl2.isValid())
    {
        Transformation goalTransformation=Transformation::identity;

        MD::Point a1=amide1.getMajorAtom()->getPosition();
        MD::Point a2=amide1.getBondSite();
        MD::Point c1=carboxyl1.getMajorAtom()->getPosition();
        MD::Point c2=carboxyl1.getBondSite();

        /* Align the two residues' average bonding site directions: */
        MD::Vector d1=Geometry::mid(amide2.getBondSite(),carboxyl2.getBondSite())-Geometry::mid(amide2.getMajorAtom()->getPosition(),carboxyl2.getMajorAtom()->getPosition());
        MD::Vector d2=Geometry::mid(a2,c2)-Geometry::mid(a1,c1);
        MD::Vector axis=Geometry::cross(d1,d2);
	dentest = (d1[0]*d1[0] + d1[1]*d1[1] + d1[2]*d1[2]) *
	      (d2[0]*d2[0] + d2[1]*d2[1] + d2[2]*d2[2]);
	if (dentest != 0.) {
	   test = (axis[0]*axis[0] + axis[1]*axis[1] +
	      axis[2]*axis[2]) / dentest;
	   if(test > .00000001) {
              MD::Scalar angle1=Math::acos(-(d1*d2)/Math::sqrt(Geometry::sqr(d1)*Geometry::sqr(d2)));
              Transformation trans1=Transformation::rotate(Transformation::Rotation::rotateAxis(axis,angle1));
              goalTransformation.leftMultiply(trans1);

              a1=trans1.transform(a1);
              a2=trans1.transform(a2);
              c1=trans1.transform(c1);
              c2=trans1.transform(c2);
              d2=trans1.transform(d2);
	      }
	   }

        /* Align the two planes formed by the residues' bonding sites: */
        MD::Vector p1=Geometry::cross(d1,carboxyl2.getBondSite()-amide2.getBondSite());
	MD::Vector plane = normalize(p1);
	MD::Point center = Geometry::mid(carboxyl2.getBondSite(), amide2.getBondSite());
	if (planepr)
	   printf("zipAntiParallel center %f %f %f\n", center[0], center[1], center[2]);
	for(int i = 0; i < 3; ++i) {
	   pl[i] = - plane[i];
	   cn[i] = center[i];
	   bb0[i] = carboxyl2.getBondSite()[i];
	   bb1[i] = amide2.getBondSite()[i];
	   }
	MD::Vector e2=c2-a2;
	MD::Vector p2=Geometry::cross(d2,e2);
        dentest = (d2[0]*d2[0] + d2[1]*d2[1] + d2[2]*d2[2]) *
           (e2[0]*e2[0] + e2[1]*e2[1] + e2[2]*e2[2]);
	if(dentest != 0.) {
	   test = (p2[0]*p2[0] + p2[1]*p2[1] + p2[2]*p2[2]) / dentest;
	   if(test > .00000001) {
              MD::Scalar angle2=Math::acos((p1*p2)/Math::sqrt(Geometry::sqr(p1)*Geometry::sqr(p2)));
              if(Geometry::cross(p2,p1)*d2<MD::Scalar(0))
                  angle2=-angle2;
              Transformation trans2=Transformation::rotate(Transformation::Rotation::rotateAxis(d2,angle2));
              goalTransformation.leftMultiply(trans2);

              a2=trans2.transform(a2);
              c2=trans2.transform(c2);
	      }
	   }

        /* Align the midpoints of the two residues' bonding sites: */
        MD::Vector dist=Geometry::mid(amide2.getBondSite(),carboxyl2.getBondSite())-Geometry::mid(a2,c2);
        Transformation trans3=Transformation::translate(dist);
        goalTransformation.leftMultiply(trans3);

                try
                        {
                        applyTransformation(goalTransformation);
                        }
                catch(DenseMatrix::RankDeficientError err)
                        {
                        /* don't do anything... */
                        fprintf(stderr,"Caught rank-deficient matrix\n");
                        }

	if(checkBondSites) {
	   e2 = amide1.getBondSite() - carboxyl2.getBondSite();
	   p2 = carboxyl1.getBondSite() - amide2.getBondSite();
	   double sum = 0;
	   for (int i = 0; i < 3; ++i)
	      sum += e2[i]*e2[i] + p2[i]*p2[i];
	   if(sum > 1.1 && lastIKResidual < .01) {
	      printf("alignment error in zipAntiParallel is %g residual %g\n",
	         sum, lastIKResidual);
	      lastIKResidual += 1.;
	      }
	   }
        return(1);
    }
		else {
		   if(prolpr) printf("zipAntiParallel %s %d %d %d %s %d %d %d is not valid\n",
		      manipulatedResidue->getPdbResidueName(),
		      manipulatedResidue->getPdbResidueIndex(), amide1.isValid(), carboxyl1.isValid(),
		      anchorResidue->getPdbResidueName(),
		      anchorResidue->getPdbResidueIndex(), amide2.isValid(), carboxyl2.isValid());
		   return(0);
		   }
    }

int zipParallel(MD::Protein::Residue* manipulatedResidue,MD::Protein::Residue* anchorResidue, float pl[3], float cn[3])
{
    typedef ProteinInteractor::Transformation Transformation;
    float test, dentest;

    /* Zip the two selected residues together: */
    MD::Protein::Dipole amide1=manipulatedResidue->getAmide();
    MD::Protein::Dipole carboxyl1=manipulatedResidue->getCarboxyl();
    MD::Protein::Residue* anchor1=anchorResidue->getPred();
    MD::Protein::Residue* anchor2=anchorResidue->getSucc();
    if(anchor1!=0&&anchor2!=0)
    {
        MD::Protein::Dipole amide2=anchor2->getAmide();
        MD::Protein::Dipole carboxyl2=anchor1->getCarboxyl();
        if(amide1.isValid()&&carboxyl1.isValid()&&amide2.isValid()&&carboxyl2.isValid())
        {
            Transformation goalTransformation=Transformation::identity;

            MD::Point a1=amide1.getBondSite();
            MD::Point c1=carboxyl1.getBondSite();
            MD::Point p1=Geometry::mid(a1,c1);
            MD::Vector d1=c1-a1;
            MD::Vector s1=(a1-amide1.getMajorAtom()->getPosition())+(c1-carboxyl1.getMajorAtom()->getPosition());
            MD::Point a2=amide2.getBondSite();
            MD::Point c2=carboxyl2.getBondSite();
            MD::Point p2=Geometry::mid(a2,c2);
            MD::Vector d2=a2-c2;
            MD::Vector s2=(a2-amide2.getMajorAtom()->getPosition())+(c2-carboxyl2.getMajorAtom()->getPosition());

            /* Align the two residues' average bonding site directions: */
            MD::Vector axis=Geometry::cross(d1,d2);
	    dentest = (d1[0]*d1[0] + d1[1]*d1[1] + d1[2]*d1[2]) *
	          (d2[0]*d2[0] + d2[1]*d2[1] + d2[2]*d2[2]);
	    if (dentest != 0.) {
	       test = (axis[0]*axis[0] + axis[1]*axis[1] +
	          axis[2]*axis[2]) / dentest;
	       if(test > .00000001) {
                  MD::Scalar angle1=Math::acos((d1*d2)/Math::sqrt(Geometry::sqr(d1)*Geometry::sqr(d2)));
                  Transformation trans1=Transformation::rotate(Transformation::Rotation::rotateAxis(axis,angle1));
                  goalTransformation.leftMultiply(trans1);

                  p1=trans1.transform(p1);
                  d1=trans1.transform(d1);
                  s1=trans1.transform(s1);
		  }
	       }

            /* Align the two planes formed by the residues' bonding sites: */
            MD::Vector pl1=Geometry::cross(d1,s1);
            MD::Vector pl2=Geometry::cross(s2,d2);
            MD::Vector axis2=Geometry::cross(pl1,pl2);
	    MD::Vector plane = normalize(pl2);
	    MD::Point center = p2;
	    if (planepr)
	       printf("zipParallel center %f %f %f\n", center[0], center[1], center[2]);
	    for(int i = 0; i < 3; ++i) {
	       pl[i] = -plane[i];
	       // cn[i] = -center[i];
	       cn[i] = center[i];
	       bb0[i] = a2[i];
	       bb1[i] = c2[i];
	       }
            dentest = (pl1[0]*pl1[0] + pl1[1]*pl1[1] + pl1[2]*pl1[2]) *
                   (pl2[0]*pl2[0] + pl2[1]*pl2[1] + pl2[2]*pl2[2]);
	    if (dentest != 0.) {
	       test = (axis2[0]*axis2[0] + axis2[1]*axis2[1] +
	          axis2[2]*axis2[2]) / dentest;
	       if(test > .00000001) {
                  MD::Scalar angle2=Math::acos((pl1*pl2)/Math::sqrt(Geometry::sqr(pl1)*Geometry::sqr(pl2)));
                  Transformation trans2=Transformation::rotate(Transformation::Rotation::rotateAxis(axis2,angle2));
                  goalTransformation.leftMultiply(trans2);

                  p1=trans2.transform(p1);
		  }
	       }

            /* Align the midpoints of the two residues' bonding sites: */
            Transformation trans3=Transformation::translate(p2-p1);
            goalTransformation.leftMultiply(trans3);

                try
                        {
                        applyTransformation(goalTransformation);
                        }
                catch(DenseMatrix::RankDeficientError err)
                        {
                        /* don't do anything... */
                        fprintf(stderr,"Caught rank-deficient matrix\n");
                        }

	    if(checkBondSites) {
	       MD::Vector e2 = amide1.getBondSite() - carboxyl2.getBondSite();
	       MD::Vector p2 = carboxyl1.getBondSite() - amide2.getBondSite();
	       double sum = 0;
	       for (int i = 0; i < 3; ++i)
	          sum += e2[i]*e2[i] + p2[i]*p2[i];
	       if(sum > 1.1 && lastIKResidual < .01) {
	          lastIKResidual += 1.;
	          printf("alignment error in zipParallel is %g\n", sum);
		  }
	       }
	    return(1);
        }
		else {
		   if(prolpr) printf("zipParallel %s %d %d %d %s %d %d %d is not valid\n",
		      manipulatedResidue->getPdbResidueName(),
		      manipulatedResidue->getPdbResidueIndex(), amide1.isValid(), carboxyl1.isValid(),
		      anchorResidue->getPdbResidueName(),
		      anchorResidue->getPdbResidueIndex(), amide2.isValid(), carboxyl2.isValid());
		   return(0);
		   }
                }
    else return(0);
    }

int altZipAntiParallel(MD::Protein::Residue* manipulatedResidue,MD::Protein::Residue* anchorResidue, float pl[3], float cn[3])
	{
	typedef ProteinInteractor::Transformation Transformation;
	float test, dentest;

	/* Zip the two selected residues together: */
	MD::Protein::Residue* anchor1=anchorResidue->getPred();
	MD::Protein::Residue* anchor2=anchorResidue->getSucc();
	MD::Protein::Residue* manip1=manipulatedResidue->getSucc();
	MD::Protein::Residue* manip2=manipulatedResidue->getPred();
	if(anchor1 == 0 || anchor2 == 0 || manip1 == 0 || manip2 == 0) return(0);
	MD::Protein::Dipole amide1=manip1->getAmide();
	MD::Protein::Dipole carboxyl1=manip2->getCarboxyl();
	MD::Protein::Dipole amide2=anchor2->getAmide();
	MD::Protein::Dipole carboxyl2=anchor1->getCarboxyl();
	if(amide1.isValid()&&carboxyl1.isValid()&&amide2.isValid()&&carboxyl2.isValid())
		{
		Transformation goalTransformation=Transformation::identity;

		MD::Point a1=amide1.getMajorAtom()->getPosition();
		MD::Point a2=amide1.getBondSite();
		MD::Point c1=carboxyl1.getMajorAtom()->getPosition();
		MD::Point c2=carboxyl1.getBondSite();

		/* Align the two residues' average bonding site directions: */
		MD::Vector d1=Geometry::mid(amide2.getBondSite(),carboxyl2.getBondSite())-Geometry::mid(amide2.getMajorAtom()->getPosition(),carboxyl2.getMajorAtom()->getPosition());
		MD::Vector d2=Geometry::mid(a2,c2)-Geometry::mid(a1,c1);
		MD::Vector axis=Geometry::cross(d1,d2);
	        dentest = (d1[0]*d1[0] + d1[1]*d1[1] + d1[2]*d1[2]) *
	              (d2[0]*d2[0] + d2[1]*d2[1] + d2[2]*d2[2]);
		if (dentest != 0.) {
		   test = (axis[0]*axis[0] + axis[1]*axis[1] +
		      axis[2]*axis[2]) / dentest;
		   if(test > .00000001) {
		      MD::Scalar angle1=Math::acos(-(d1*d2)/Math::sqrt(Geometry::sqr(d1)*Geometry::sqr(d2)));
		      Transformation trans1=Transformation::rotate(Transformation::Rotation::rotateAxis(axis,angle1));
		      goalTransformation.leftMultiply(trans1);

		      a1=trans1.transform(a1);
		      a2=trans1.transform(a2);
		      c1=trans1.transform(c1);
		      c2=trans1.transform(c2);
		      d2=trans1.transform(d2);
		      }
		   }

		/* Align the two planes formed by the residues' bonding sites: */
		MD::Vector e1 = carboxyl2.getBondSite()-amide2.getBondSite();
		MD::Vector p1=Geometry::cross(d1,carboxyl2.getBondSite()-amide2.getBondSite());
	        MD::Vector plane = normalize(p1);
	        MD::Point center = Geometry::mid(carboxyl2.getBondSite(), amide2.getBondSite());
	        if (planepr)
	           printf("altZipAntiParallel center %f %f %f\n", center[0], center[1], center[2]);
		for(int i = 0; i < 3; ++i) {
		   pl[i] = plane[i];
		   cn[i] = center[i];
	           bb0[i] = carboxyl2.getBondSite()[i];
	           bb1[i] = amide2.getBondSite()[i];
		   }
		MD::Vector e2=c2-a2;
		MD::Vector p2=Geometry::cross(d2,e2);
	        dentest = (d2[0]*d2[0] + d2[1]*d2[1] + d2[2]*d2[2]) *
	           (e2[0]*e2[0] + e2[1]*e2[1] + e2[2]*e2[2]);
		if(dentest != 0.) {
		   test = (p2[0]*p2[0] + p2[1]*p2[1] +
		      p2[2]*p2[2]) / dentest;
		   if(dentest > .00000001) {
		      MD::Scalar angle2=Math::acos((p1*p2)/Math::sqrt(Geometry::sqr(p1)*Geometry::sqr(p2)));
		      if(Geometry::cross(p2,p1)*d2<MD::Scalar(0))
			      angle2=-angle2;
		      Transformation trans2=Transformation::rotate(Transformation::Rotation::rotateAxis(d2,angle2));
		      goalTransformation.leftMultiply(trans2);

		      a2=trans2.transform(a2);
		      c2=trans2.transform(c2);
		      }
		   }

		/* Align the midpoints of the two residues' bonding sites: */
		MD::Vector dist=Geometry::mid(amide2.getBondSite(),carboxyl2.getBondSite())-Geometry::mid(a2,c2);
		Transformation trans3=Transformation::translate(dist);
		goalTransformation.leftMultiply(trans3);

                try
                        {
                        applyTransformation(goalTransformation);
                        }
                catch(DenseMatrix::RankDeficientError err)
                        {
                        /* don't do anything... */
                        fprintf(stderr,"Caught rank-deficient matrix\n");
                        }

	       if(checkBondSites) {
	          e2 = amide1.getBondSite() - carboxyl2.getBondSite();
	          p2 = carboxyl1.getBondSite() - amide2.getBondSite();
	          double sum = 0;
	          for (int i = 0; i < 3; ++i)
	             sum += e2[i]*e2[i] + p2[i]*p2[i];
	          if(sum > 1.1 && lastIKResidual < .01) {
		     lastIKResidual += 1.;
	             printf("alignment error in altZipAntiParallel is %g\n", sum);
		     }
	          }
		return(1);
		}
		else {
		   if(prolpr) printf("altZipAntiParallel %s %d %d %d %s %d %d %d is not valid\n",
		      manipulatedResidue->getPdbResidueName(),
		      manipulatedResidue->getPdbResidueIndex(), amide1.isValid(), carboxyl1.isValid(),
		      anchorResidue->getPdbResidueName(),
		      anchorResidue->getPdbResidueIndex(), amide2.isValid(), carboxyl2.isValid());
		   return(0);
		   }
    }

int altZipParallel(MD::Protein::Residue* manipulatedResidue,MD::Protein::Residue* anchorResidue, float pl[3], float cn[3])
	{
	typedef ProteinInteractor::Transformation Transformation;
	float test, dentest;

	/* Zip the two selected residues together: */
	MD::Protein::Dipole amide2=anchorResidue->getAmide();
	MD::Protein::Dipole carboxyl2=anchorResidue->getCarboxyl();
	MD::Protein::Residue* manip1=manipulatedResidue->getPred();
	MD::Protein::Residue* manip2=manipulatedResidue->getSucc();
	if(manip1!=0&&manip2!=0)
		{
		MD::Protein::Dipole amide1=manip2->getAmide();
		MD::Protein::Dipole carboxyl1=manip1->getCarboxyl();
		if(amide1.isValid()&&carboxyl1.isValid()&&amide2.isValid()&&carboxyl2.isValid())
			{
			Transformation goalTransformation=Transformation::identity;

			MD::Point a1=amide1.getBondSite();
			MD::Point c1=carboxyl1.getBondSite();
			MD::Point p1=Geometry::mid(a1,c1);
			MD::Vector d1=c1-a1;
			MD::Vector s1=(a1-amide1.getMajorAtom()->getPosition())+(c1-carboxyl1.getMajorAtom()->getPosition());
			MD::Point a2=amide2.getBondSite();
			MD::Point c2=carboxyl2.getBondSite();
			MD::Point p2=Geometry::mid(a2,c2);
			MD::Vector d2=a2-c2;
			MD::Vector s2=(a2-amide2.getMajorAtom()->getPosition())+(c2-carboxyl2.getMajorAtom()->getPosition());

			/* Align the two residues' average bonding site directions: */
			MD::Vector axis=Geometry::cross(d1,d2);

	                dentest = (d1[0]*d1[0] + d1[1]*d1[1] + d1[2]*d1[2]) *
	                   (d2[0]*d2[0] + d2[1]*d2[1] + d2[2]*d2[2]);
			if (dentest != 0.) {
			   test = (axis[0]*axis[0] + axis[1]*axis[1] +
			      axis[2]*axis[2]) / dentest;
			   if(test > .00000001) {
			      MD::Scalar angle1=Math::acos((d1*d2)/Math::sqrt(Geometry::sqr(d1)*Geometry::sqr(d2)));
			      Transformation trans1=Transformation::rotate(Transformation::Rotation::rotateAxis(axis,angle1));
			      goalTransformation.leftMultiply(trans1);
			      p1=trans1.transform(p1);
			      d1=trans1.transform(d1);
			      s1=trans1.transform(s1);
			      }
			   }

			/* Align the two planes formed by the residues' bonding sites: */
			MD::Vector pl1=Geometry::cross(d1,s1);
			MD::Vector pl2=Geometry::cross(s2,d2);
			MD::Vector axis2=Geometry::cross(pl1,pl2);
	                MD::Vector plane = normalize(pl2);
	                MD::Point center = p2;
	        if (planepr)
	           printf("altZipParallel center %f %f %f\n", center[0], center[1], center[2]);
			for(int i = 0; i < 3; ++i) {
			   pl[i] = plane[i];
			   cn[i] = center[i];
	                   bb0[i] = a2[i];
	                   bb1[i] = c2[i];
			   }
	                dentest = (pl1[0]*pl1[0] + pl1[1]*pl1[1] + pl1[2]*pl1[2]) *
	                   (pl2[0]*pl2[0] + pl2[1]*pl2[1] + pl2[2]*pl2[2]);
			if (dentest != 0.) {
			   test = (axis2[0]*axis2[0] + axis2[1]*axis2[1] +
			      axis2[2]*axis2[2]) / dentest;
			   if(test > .00000001) {
			      MD::Scalar angle2=Math::acos((pl1*pl2)/Math::sqrt(Geometry::sqr(pl1)*Geometry::sqr(pl2)));
			      Transformation trans2=Transformation::rotate(Transformation::Rotation::rotateAxis(axis2,angle2));
			      goalTransformation.leftMultiply(trans2);

			      p1=trans2.transform(p1);
			      }
			   }

			/* Align the midpoints of the two residues' bonding sites: */
			Transformation trans3=Transformation::translate(p2-p1);
			goalTransformation.leftMultiply(trans3);

                	try
                        	{
                        	applyTransformation(goalTransformation);
                        	}
                	catch(DenseMatrix::RankDeficientError err)
                        	{
                        	/* don't do anything... */
                        	fprintf(stderr,"Caught rank-deficient matrix\n");
                        	}

	if(checkBondSites) {
	   MD::Vector e2 = amide1.getBondSite() - carboxyl2.getBondSite();
	   MD::Vector p2 = carboxyl1.getBondSite() - amide2.getBondSite();
	   double sum = 0;
	   for (int i = 0; i < 3; ++i)
	      sum += e2[i]*e2[i] + p2[i]*p2[i];
	   if(sum > 1.1 && lastIKResidual < .01) {
	      lastIKResidual += 1.;
	      printf("alignment error in altZipParallel is %g\n", sum);
	      }
	   }
                        return(1);
                }
                else {
		   if(prolpr) printf("altZipParallel %s %d %d %d %s %d %d %d is not valid\n",
		      manipulatedResidue->getPdbResidueName(),
		      manipulatedResidue->getPdbResidueIndex(), amide1.isValid(), carboxyl1.isValid(),
		      anchorResidue->getPdbResidueName(),
		      anchorResidue->getPdbResidueIndex(), amide2.isValid(), carboxyl2.isValid());
                   return(0);
                   }
                }
        else return(0);
        }

int main(int argc, char* argv[]) {
	try {
		char** appDefaults = 0;
		VRProteinShop app(argc, argv, appDefaults);
		app.run();
	} catch (std::runtime_error err) {
		std::cerr << "Caught exception " << err.what() << std::endl;
		return 1;
	}

	return 0;
}
