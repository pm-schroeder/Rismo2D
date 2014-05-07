#-------------------------------------------------
#
# Project created by QtCreator 2012-01-27T09:35:46
#
#-------------------------------------------------

QT += core
QT -= gui

TARGET = RismoQt

CONFIG += console
CONFIG -= app_bundle

TEMPLATE = app


SOURCES += \
    Velzen.cpp \
    VeloGrad.cpp \
    Update.cpp \
    Type.cpp \
    Turbulence.cpp \
    Triangle.cpp \
    Transform.cpp \
    Timeint.cpp \
    Time.cpp \
    Subdom.cpp \
    Statist.cpp \
    Square.cpp \
    Solver.cpp \
    Solve.cpp \
    Smooth.cpp \
    SlipFlow.cpp \
    Shape.cpp \
    SetEqno.cpp \
    SetBdKD.cpp \
    Sed.cpp \
    Section.cpp \
    Scale.cpp \
    Rotate.cpp \
    Rot2D.cpp \
    Report.cpp \
    ReorderElem.cpp \
    Reorder.cpp \
    Project.cpp \
    Preco_ilut.cpp \
    Preco_ilu0.cpp \
    Phi2D.cpp \
    P_fgmresd.cpp \
    P_bcgstabd.cpp \
    Node.cpp \
    Model.cpp \
    Memory.cpp \
    Main.cpp \
    Lumped.cpp \
    Locate.cpp \
    Line.cpp \
    Lindner.cpp \
    LastNode.cpp \
    Interpol.cpp \
    InitS.cpp \
    Init.cpp \
    IndexMat.cpp \
    Grid.cpp \
    Frontm.cpp \
    Front.cpp \
    Fromat.cpp \
    Friction.cpp \
    EqsUVS2D_LV.cpp \
    EqsUVS2D.cpp \
    EqsSL2D.cpp \
    EqsPPE2D.cpp \
    EqsKL2D.cpp \
    EqsKD2D.cpp \
    EqsK2D.cpp \
    EqsDz.cpp \
    EqsDisp.cpp \
    EqsD2D.cpp \
    EqsBL2D.cpp \
    Eqs.cpp \
    Elem.cpp \
    EddyDisp.cpp \
    DryRewet.cpp \
    DoFriction.cpp \
    DoDryRew.cpp \
    Dispersion.cpp \
    Cycle.cpp \
    Curv2D.cpp \
    CRSMat.cpp \
    Courant.cpp \
    ContinBSL.cpp \
    Contin.cpp \
    Connect.cpp \
    Compute.cpp \
    CoefsUVS2D_LV.cpp \
    CoefsUVS2D.cpp \
    CoefsSL2D.cpp \
    CoefsPPE2D.cpp \
    CoefsKL2D.cpp \
    CoefsKD2D.cpp \
    CoefsK2D.cpp \
    CoefsDz.cpp \
    CoefsDisp.cpp \
    CoefsD2D.cpp \
    CoefsBL2D.cpp \
    Check.cpp \
    Bound.cpp \
    Bicgstab.cpp \
    BconSet.cpp \
    BconLine.cpp \
    Bcon.cpp \
    Assemble.cpp \
    Asciifile.cpp \
    ArFact.cpp \
    CoefsUVS2D_AI.cpp \
    EqsUVS2D_AI.cpp \
    CoefsUVS2D_TM.cpp \
    EqsUVS2D_TM.cpp \
    EqsUVS2D_TMAI.cpp \
    CoefsUVS2D_TMAI.cpp

HEADERS += \
    Vars.h \
    Type.h \
    Tools.h \
    Tmserhead.h \
    Timeint.h \
    Time.h \
    Subdom.h \
    Statist.h \
    Solver.h \
    Shape.h \
    Sed.h \
    Section.h \
    Scale.h \
    Report.h \
    Reorder.h \
    Project.h \
    Precon.h \
    Preco_ilut.h \
    Preco_ilu0.h \
    P_fgmresd.h \
    P_bcgstabd.h \
    Parms.h \
    Node.h \
    Model.h \
    Memory.h \
    Grid.h \
    Frontm.h \
    Front.h \
    Fromat.h \
    EqsUVS2D_LV.h \
    EqsUVS2D.h \
    EqsSL2D.h \
    EqsPPE2D.h \
    EqsKL2D.h \
    EqsKD2D.h \
    EqsK2D.h \
    EqsDz.h \
    EqsDisp.h \
    EqsD2D.h \
    EqsBL2D.h \
    Eqs.h \
    Elem.h \
    Defs.h \
    Datkey.h \
    CRSMat.h \
    Bicgstab.h \
    Bcon.h \
    Asciifile.h \
    EqsUVS2D_AI.h \
    EqsUVS2D_TM.h \
    EqsUVS2D_TMAI.h

#INCLUDEPATH += /usr/include/mpi

#unix:!macx:!symbian: LIBS += -lmpi -lmpi++
