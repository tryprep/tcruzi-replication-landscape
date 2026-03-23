#!/usr/bin/env bash

###################################################################################################
#                                        BCFtoolsDriver                                           #
#                                                                                                 #
#  Description:                                                                                   #
#  This script orchestrates the execution of all pre-processing, processing, and post-processing  #
#  steps in the comprehensive SNP analysis pipeline, utilizing variant calls from FreeBayes.      #
#                                                                                                 #
#  As the master controller, it directs the investigation of SNP density and distribution across  #
#  DNA replication origins and termini, multigene families, and specific chromosomes. It also     #
#  triggers advanced downstream analyses, including SnpEff functional annotation, dN/dS ratio     #
#  calculations, and Loss of Heterozygosity (LOH) evaluations, using a combination of BCFtools,   #
#  custom Python, and C++ scripts.                                                                #
#                                                                                                 #
#  Usage:                                                                                         #
#  saveCommand script/BCFtoolsDriver.bash                                                         #
#                                                                                                 #
#  Copyleft (ɔ) 2025 by Thiago Andrade Franco                                                     #
#                       <thiago.franco@fundacaobutantan.org.br>                                   #
#                                                                                                 #
#  23/04/2025: First version.                                                                     #
###################################################################################################

# Loading environment variables about this experiment.
. script/BCFtoolsVariables.bash

# Running genomeCompartmentComposition pre-processing, main and post-processing scripts.
saveCommand script/BCFtoolsPre.bash
saveCommand script/BCFtoolsMain.bash
saveCommand script/BCFtoolsMain2.bash
saveCommand script/BCFtoolsMain3.bash
saveCommand script/BCFtoolsMain4.bash

exit 0
