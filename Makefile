OBJ  = sources/ArFact.o             sources/Asciifile.o\
       sources/Assemble.o           sources/Bcon.o\
       sources/BconLine.o           sources/BconSet.o\
       sources/Bicgstab.o           sources/Bound.o\
       sources/Check.o              sources/CoefsBL2D.o\
       sources/CoefsD2D.o           sources/CoefsDisp.o\
       sources/CoefsDz.o            sources/CoefsK2D.o\
       sources/CoefsKD2D.o          sources/CoefsKL2D.o\
       sources/CoefsPPE2D.o         sources/CoefsSL2D.o\
       sources/CoefsUVS2D.o         sources/CoefsUVS2D_AI.o\
       sources/CoefsUVS2D_TM.o      sources/CoefsUVS2D_TMAI.o\
       sources/CoefsUVS2D_ME.o      sources/CoefsUVS2D_ME_AI.o\
       sources/CoefsUVS2D_ME_TM.o   sources/CoefsUVS2D_ME_TMAI.o\
       sources/Compute.o            sources/Connect.o\
       sources/ContinBSL.o          sources/Contin.o\
       sources/Courant.o            sources/CRSMat.o\
       sources/Curv2D.o             sources/Cycle.o\
       sources/Dispersion.o         sources/DoDryRew.o\
       sources/DoFriction.o         sources/DryRewet.o\
       sources/EddyDisp.o           sources/Elem.o\
       sources/EqsBL2D.o            sources/Eqs.o\
       sources/EqsD2D.o             sources/EqsDisp.o\
       sources/EqsDz.o              sources/EqsK2D.o\
       sources/EqsKD2D.o            sources/EqsKL2D.o\
       sources/EqsPPE2D.o           sources/EqsSL2D.o\
       sources/EqsUVS2D.o           sources/EqsUVS2D_AI.o\
       sources/EqsUVS2D_TM.o        sources/EqsUVS2D_TMAI.o\
       sources/EqsUVS2D_ME.o        sources/EqsUVS2D_ME_AI.o\
       sources/EqsUVS2D_ME_TM.o     sources/EqsUVS2D_ME_TMAI.o\
       sources/Friction.o           sources/Fromat.o\
       sources/Front.o              sources/Frontm.o\
       sources/Grid.o\
       sources/IndexMat.o           sources/Init.o\
       sources/InitS.o              sources/Interpol.o\
       sources/LastNode.o           sources/Lindner.o\
       sources/Line.o               sources/Locate.o\
       sources/Lumped.o             sources/Main.o\
       sources/Memory.o             sources/Model.o\
       sources/Node.o               sources/P_bcgstabd.o\
       sources/P_fgmresd.o          sources/Phi2D.o\
       sources/Preco_ilu0.o         sources/Preco_ilut.o\
       sources/Project.o            sources/Reorder.o\
       sources/ReorderElem.o        sources/Report.o\
       sources/Rot2D.o              sources/Rotate.o\
       sources/Scale.o              sources/Section.o\
       sources/Sed.o                sources/SetBdKD.o\
       sources/SetEqno.o            sources/Shape.o\
       sources/SlipFlow.o           sources/Smooth.o\
       sources/Solve.o              sources/Solver.o\
       sources/Square.o             sources/Statist.o\
       sources/Subdom.o             sources/Time.o\
       sources/Timeint.o            sources/Transform.o\
       sources/Triangle.o           sources/Turbulence.o\
       sources/Type.o               sources/Update.o\
       sources/VeloGrad.o           sources/Velzen.o

# -----------------------------------------------------------------------------------------
# Namen

PROG = rismo_40600

# ---------------------------------------------------------
# Optionen

COMP = g++

COPT =
LOPT =

# ---------------------------------------------------------

.SUFFIXES : .o .cpp .cpp~
.cpp.o:

# ---------------------------------------------------------

$(PROG) : $(OBJ)
	$(COMP) $(OBJ) $(LOPT) -o $(PROG)

$(OBJ) :
	$(COMP) $(COPT) -c $*.cpp -o $*.o
