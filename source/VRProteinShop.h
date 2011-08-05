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

#ifndef VRPROTEINSHOP_INCLUDED
#define VRPROTEINSHOP_INCLUDED

#include <vector>
#include <Misc/ConfigurationFile.h>
#include <Threads/Thread.h>
#include <Threads/Mutex.h>
#include <Threads/Cond.h>
#include <Geometry/Ray.h>
#include <GL/gl.h>
#include <GL/GLObject.h>
#include <GLMotif/RadioBox.h>
#include <GLMotif/ToggleButton.h>
#include <Vrui/Tools/LocatorTool.h>
#include <Vrui/LocatorToolAdapter.h>
#include <Vrui/Tools/DraggingTool.h>
#include <Vrui/DraggingToolAdapter.h>
#include <Vrui/ToolManager.h>
#include <Vrui/Application.h>

#include <Protein/PaletteRenderer.h>

/* Forward declarations: */
namespace Comm {
class MulticastPipe;
}
class GLColorMap;
namespace GLMotif {
class Popup;
class PopupMenu;
class PopupWindow;
}
namespace MD {
class Protein;
class ProteinRenderer;
}
class UndoBuffer;
class ProteinInteractor;
class EnergyCalculator;
class ClientServerPipe;

struct VRProteinShop:public Vrui::Application,public GLObject
	{
	/* Embedded classes: */
	private:
	enum SelectionMode
		{
		SELECT_SECONDARYSTRUCTURE,SELECT_SINGLERESIDUE,SELECT_TOGGLEACTIVECOIL
		};
	
	enum AutoBondingMode
		{
		AUTOBOND_NONE,AUTOBOND_PARALLEL,AUTOBOND_ANTIPARALLEL
		};
	
	enum UpdateDirection
		{
		UPDATE_LEFT,UPDATE_RIGHT
		};
	
	class Locator:public Vrui::LocatorToolAdapter // Base class for locator tool adapters
		{
		/* Elements: */
		protected:
		VRProteinShop* application;
		
		/* Constructors and destructors: */
		public:
		Locator(Vrui::LocatorTool* sTool,VRProteinShop* sApplication)
			:Vrui::LocatorToolAdapter(sTool),
			 application(sApplication)
			{
			};
		};
	
	typedef std::vector<Locator*> LocatorList;
	
	class StructureLocator:public Locator
		{
		/* Constructors and destructors: */
		public:
		StructureLocator(Vrui::LocatorTool* sTool,VRProteinShop* sApplication)
			:Locator(sTool,sApplication)
			{
			};
		
		/* Methods: */
		virtual void buttonPressCallback(Vrui::LocatorTool::ButtonPressCallbackData* cbData);
		};
	
	class Dragger:public Vrui::DraggingToolAdapter // Base class for dragging tool adapters
		{
		/* Elements: */
		protected:
		VRProteinShop* application;
		
		/* Constructors and destructors: */
		public:
		Dragger(Vrui::DraggingTool* sTool,VRProteinShop* sApplication)
			:Vrui::DraggingToolAdapter(sTool),
			 application(sApplication)
			{
			};
		};
	
	class StructureDragger:public Dragger
		{
		/* Elements: */
		private:
		bool isDragging; // Flag if the dragger is active
		bool rayBased; // Flag if the dragging is to use the ray dragging method
		Vrui::Ray initialRay; // Equation of initial dragging ray
		
		/* Constructors and destructors: */
		public:
		StructureDragger(Vrui::DraggingTool* sTool,VRProteinShop* sApplication)
			:Dragger(sTool,sApplication),
			 isDragging(false)
			{
			};
		
		/* Methods: */
		virtual void dragStartCallback(Vrui::DraggingTool::DragStartCallbackData* cbData);
		virtual void dragCallback(Vrui::DraggingTool::DragCallbackData* cbData);
		virtual void dragEndCallback(Vrui::DraggingTool::DragEndCallbackData* cbData);
		};
	
	struct DataItem:public GLObject::DataItem
		{
		/* Elements: */
		public:
		GLuint graphicsDisplayListId;
		unsigned int graphicsVersion;
		
		/* Constructors and destructors: */
		DataItem(void)
			:graphicsDisplayListId(glGenLists(1)),graphicsVersion(0)
			{
			};
		~DataItem(void)
			{
			glDeleteLists(graphicsDisplayListId,1);
			};
		};
	
	typedef std::vector<Dragger*> DraggerList;
	
	friend class StructureLocator;
	friend class StructureDragger;
	
	/* Elements: */
	private:

	bool validPdbForEnergy;

	/* Protein stuff: */
	mutable Threads::Mutex proteinMutex; // Mutex protecting the protein
	MD::Protein* protein; // The current protein
	unsigned int proteinID;
	unsigned int graphicsVersion; // "Graphics" version of the protein; incremented every time the protein is changed
	EnergyCalculator* energyCalculator;
	
	/* Rendering stuff: */
	MD::ProteinRenderer* proteinRenderer; // Renderer attached to the current protein
	bool drawDensity; // Flag whether to render the density distribution
	PaletteRenderer* densityRenderer; // A volume renderer to render density distributions
	GLColorMap* densityPalette; // Color map to render density distributions
	VolumeRenderer::Vector viewDirection; // Current view direction
	
	/* Interaction state: */
	SelectionMode selectionMode; // Current selection mode
	AutoBondingMode autoBondingMode; // Current auto-bonding mode
	UpdateDirection updateDirection; // Current update direction
	UndoBuffer * undoBuffer;
	ProteinInteractor* interactor; // Pointer to protein interaction object
	Comm::MulticastPipe* interactorPipe; // Multicast pipe to synchronize IK updates across a distributed rendering cluster
	Threads::Cond ikUpdateRequestedCond;
	volatile bool ikUpdateRequested;
	Threads::Cond ikUpdateDoneCond;
	volatile bool ikUpdateDone;
	volatile bool ikUpdatePosted;
	Threads::Thread ikUpdateThread;
	
	/* Vrui stuff: */
	LocatorList locators; // List of active locator objects
	DraggerList draggers; // List of active dragger objects
	GLMotif::PopupMenu* mainMenu; // The program's main menu
	GLMotif::ToggleButton* globalDrawAtomsToggle;
	GLMotif::ToggleButton* globalDrawBondsToggle;
	GLMotif::ToggleButton* globalDrawBackboneToggle;
	GLMotif::ToggleButton* globalDrawCartoonToggle;
	GLMotif::ToggleButton* globalDrawHydrogenBondsToggle;
	GLMotif::ToggleButton* globalDrawHydrogenBondSitesToggle;
	GLMotif::ToggleButton* globalDrawHydrogenCagesToggle;
	GLMotif::ToggleButton* globalDrawAtomCollisionsToggle;
	GLMotif::ToggleButton* globalDrawDensityToggle;
	GLMotif::PopupWindow* interactionDialog; // The interaction dialog
	GLMotif::RadioBox* autoBondingBox;
	
	/* Private methods: */
	GLMotif::Popup* createGlobalDisplayMenu(void);
	GLMotif::PopupMenu* createMainMenu(void);
	GLMotif::PopupWindow* createInteractionDialog(void);
	void updateMenuState(void);
	void* ikUpdateThreadMethod(void);
	
	/* Constructors and destructors: */
	public:
	VRProteinShop(int& argc,char**& argv,char**& appDefaults);
	virtual ~VRProteinShop(void);
	
	/* Methods: */
	void addProtein (MD::Protein* newProtein, int Id);
	virtual void initContext(GLContextData& contextData) const;
	void loadProteins(const char *filename);
	MD::Protein * loadPdb (const char *filename);
	virtual void toolCreationCallback(Vrui::ToolManager::ToolCreationCallbackData* cbData);
	virtual void toolDestructionCallback(Vrui::ToolManager::ToolDestructionCallbackData* cbData);
	virtual void frame(void);
	virtual void display(GLContextData& contextData) const;
	void centerDisplayCallback(Misc::CallbackData* cbData);
	void menuToggleSelectCallback(GLMotif::ToggleButton::ValueChangedCallbackData* cbData);
	void showInteractionDialogValueChangedCallback(GLMotif::ToggleButton::ValueChangedCallbackData* cbData);
	void selectionModeValueChangedCallback(GLMotif::RadioBox::ValueChangedCallbackData* cbData);
	void autoBondingValueChangedCallback(GLMotif::RadioBox::ValueChangedCallbackData* cbData);
	void resetActiveCoilRegionsCallback(Misc::CallbackData* cbData);
	void ikUpdateDirectionValueChangedCallback(GLMotif::RadioBox::ValueChangedCallbackData* cbData);
	};

#endif
