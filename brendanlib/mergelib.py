import numpy as np
import glob,subprocess,os,sys
import haloutils as htils

def make_music_submission_script(runpath,cfgname,job_name,runsingle):
    f = open(runpath + "/smusic",'w')
    f.write("#!/bin/bash \n")
    f.write("#SBATCH -o music.o%j \n")
    f.write("#SBATCH -e music.e%j \n")
    f.write("#SBATCH -N 1\n")
    f.write("#SBATCH --exclusive\n")
    
    f.write('#SBATCH -p AMD64\n')
    ncores = 64
    if runsingle:
        f.write("#SBATCH -J " + job_name + "s\n")
    else:
        f.write("#SBATCH -J " + job_name + "\n")

    f.write("\n")
    f.write("export OMP_NUM_THREADS=" + str(ncores) + "\n")
    f.write("\n")
    f.write("cd " + runpath + "\n")
    f.write("source ~/.bashrc\n")
    f.write("\n")
    if runsingle:
        f.write("./MUSIC_singleprec ./" + cfgname + ".conf 1>OUTPUTmusic 2>ERRORmusic\n")
    else:
        f.write("./MUSIC ./" + cfgname + ".conf 1>OUTPUTmusic 2>ERRORmusic\n")
    f.close()

def make_arepo_submission_script(runpath,job_name,restart_flag,run_amd):
    f = open(runpath + "/sarepo",'w')
    f.write('#!/bin/bash\n')
    f.write('#SBATCH -o arepo.o%j\n')
    f.write('#SBATCH -e arepo.e%j\n')
    f.write('#SBATCH -J '+ job_name + '\n')
    
    if run_amd:
        queue = "AMD64"
        ncores = 64
        core_setup = '-n ' + str(ncores)
    else:
        queue = "HyperNodes"
        ncores = 120
        core_setup = '-n ' + str(ncores)

    f.write('#SBATCH -p ' + queue + '\n')
    f.write('#SBATCH ' + core_setup + '\n')

    f.write("\n")
    f.write("cd " + runpath + "\n")
    f.write("\n")

    if not restart_flag:
            if queue != "HyperNodes":
                f.write('mpirun -np ' + str(ncores) +  ' ./Arepo param.txt 1>OUTPUT 2>ERROR\n')
            else:
                f.write('mpirun -np ' + str(ncores) +  ' ./Arepo param.txt 1>OUTPUT 2>ERROR\n')
    else:
            if queue != "HyperNodes":
                f.write('mpirun -np ' + str(ncores) +  ' ./Arepo param.txt 1 1>>OUTPUT 2>ERROR\n')
            else:
                f.write('mpirun -np ' + str(ncores) +  ' ./Arepo param.txt 1 1>>OUTPUT 2>ERROR\n')
    f.close()

def make_arepo_parameter_file_dm(runpath,job_name,run_amd,arepofilespath):
    f = open(runpath + "/param.txt",'w')

    f.write("%----  Relevant files \n")
    f.write("InitCondFile        ./ics\n")
    f.write("OutputDir           ./outputs/\n")
    f.write("SnapshotFileBase    snap\n")
    f.write("OutputListFilename  ExpansionList\n")
    f.write("\n")
    f.write("%---- File formats\n")
    f.write("ICFormat           3 \n")
    f.write("SnapFormat         3 \n")
    f.write("\n")
    f.write("%---- CPU-time limits\n")
    f.write("TimeLimitCPU              2592000   % in seconds\n")
    f.write("CpuTimeBetRestartFile      10800    % in seconds\n")
    f.write("ResubmitOn        0\n")
    f.write("ResubmitCommand   my-scriptfile \n")
    f.write("\n")
    f.write("%----- Memory alloction\n")
    if run_amd:
        f.write("MaxMemSize        3800 \n")
    else:
        f.write("MaxMemSize        1000 \n")
    f.write("\n")
    f.write("%---- Caracteristics of run\n")
    f.write("TimeBegin           0.0078125      % Begin of the simulation\n")
    f.write("TimeMax             1.0            % End of the simulation\n")
    f.write("\n")
    f.write("%---- Basic code options that set the type of simulation\n")
    f.write("ComovingIntegrationOn    1 \n")
    f.write("PeriodicBoundariesOn     1\n")
    f.write("CoolingOn        0\n")
    f.write("StarformationOn      0 \n")
    f.write("\n")
    f.write("%---- Cosmological parameters\n")
    f.write("Omega0                0.3175\n")
    f.write("OmegaLambda           0.6825\n")
    f.write("OmegaBaryon           0.0\n")
    f.write("HubbleParam           0.6711\n")
    f.write("BoxSize               100000.0\n")
    f.write("\n")
    f.write("%---- Output frequency and output paramaters\n")
    f.write("OutputListOn              1 \n")
    f.write("TimeBetSnapshot           0.0\n")
    f.write("TimeOfFirstSnapshot       0.0\n")
    f.write("TimeBetStatistics         0.01\n")
    f.write("NumFilesPerSnapshot       8 \n")
    f.write("NumFilesWrittenInParallel 8\n")
    f.write("\n")
    f.write("%---- Accuracy of time integration\n")
    f.write("TypeOfTimestepCriterion  0 \n")
    f.write("ErrTolIntAccuracy        0.012 \n")
    f.write("CourantFac               0.3 \n")
    f.write("MaxSizeTimestep          0.0025\n")
    f.write("MinSizeTimestep          0.0 \n")
    f.write("\n")
    f.write("%---- Treatment of empty space and temperature limits\n")
    f.write("InitGasTemp            350.0\n")
    f.write("MinGasTemp             15.0\n")
    f.write("MinimumDensityOnStartUp               1.0e-20\n")
    f.write("LimitUBelowThisDensity                0.0\n")
    f.write("LimitUBelowCertainDensityToThisValue  0.0\n")
    f.write("MinEgySpec             0\n")
    f.write("\n")
    f.write("%---- Tree algorithm, force accuracy, domain update frequency\n")
    f.write("TypeOfOpeningCriterion            1\n")
    f.write("ErrTolTheta                       0.7\n")
    f.write("ErrTolForceAcc                    0.0025\n")
    f.write("MultipleDomains                   8 \n")
    f.write("TopNodeFactor                     5\n")
    f.write("ActivePartFracForNewDomainDecomp  0.01\n")
    f.write("\n")
    f.write("%---- Initial density estimate\n")
    f.write("DesNumNgb              64\n")
    f.write("MaxNumNgbDeviation     1 \n")
    f.write("\n")
    f.write("%---- System of units\n")
    f.write("UnitLength_in_cm         3.085678e21        ;  1.0 Kpc \n")
    f.write("UnitMass_in_g            1.989e43           ;  1.0e10 solar masses\n")
    f.write("UnitVelocity_in_cm_per_s 1e5                ;  1 km/sec\n")
    f.write("GravityConstantInternal  0\n")
    f.write("\n")
    f.write("%---- Gravitational softening lengths\n")
    f.write("GasSoftFactor      1.5\n")

    if "LX11" in runpath:
        f.write("SofteningComovingType0    0.0\n")
        f.write("SofteningComovingType1    0.610352\n")
        f.write("SofteningComovingType2    2.441406\n")
        f.write("SofteningComovingType3    4.882813\n")
        f.write("SofteningComovingType4    195.313\n")
        f.write("SofteningComovingType5    390.625\n")

        f.write("SofteningMaxPhysType0     0.0\n")
        f.write("SofteningMaxPhysType1     0.610352\n")
        f.write("SofteningMaxPhysType2     2.441406\n")
        f.write("SofteningMaxPhysType3     4.882813\n")
        f.write("SofteningMaxPhysType4     19.5313\n")
        f.write("SofteningMaxPhysType5     39.0625\n")

    if "LX12" in runpath:
        f.write("SofteningComovingType0    0.0\n")
        f.write("SofteningComovingType1    0.305175\n")
        f.write("SofteningComovingType2    1.220703\n")
        f.write("SofteningComovingType3    2.441406\n")
        f.write("SofteningComovingType4    9.76563\n")
        f.write("SofteningComovingType5    19.5313\n")

        f.write("SofteningMaxPhysType0     0.0\n")
        f.write("SofteningMaxPhysType1     0.305175\n")
        f.write("SofteningMaxPhysType2     1.220703\n")
        f.write("SofteningMaxPhysType3     2.441406\n")
        f.write("SofteningMaxPhysType4     9.76563\n")
        f.write("SofteningMaxPhysType5     19.5313\n")

    if "LX13" in runpath:
        f.write("SofteningComovingType0    0.0\n")
        f.write("SofteningComovingType1    0.15258789\n")
        f.write("SofteningComovingType2    0.61035156\n")
        f.write("SofteningComovingType3    1.22070313\n")
        f.write("SofteningComovingType4    4.88281\n")
        f.write("SofteningComovingType5    9.76563\n")
        
        f.write("SofteningMaxPhysType0     0.0\n")
        f.write("SofteningMaxPhysType1     0.15258789\n")
        f.write("SofteningMaxPhysType2     0.61035156\n")
        f.write("SofteningMaxPhysType3     1.22070313\n")
        f.write("SofteningMaxPhysType4     4.88281\n")
        f.write("SofteningMaxPhysType5     9.76563\n")

    if "LX14" in runpath:
        f.write("SofteningComovingType0    0.0762939453125\n")
        f.write("SofteningComovingType1    0.0762939453125\n")
        f.write("SofteningComovingType2    0.30517578\n")
        f.write("SofteningComovingType3    0.61035156\n")
        f.write("SofteningComovingType4    2.44141\n")
        f.write("SofteningComovingType5    4.88281\n")
        
        f.write("SofteningMaxPhysType0     0.0762939453125\n")
        f.write("SofteningMaxPhysType1     0.0762939453125\n")
        f.write("SofteningMaxPhysType2     0.30517578\n")
        f.write("SofteningMaxPhysType3     0.61035156\n")
        f.write("SofteningMaxPhysType4     2.44141\n")
        f.write("SofteningMaxPhysType5     4.88281\n")

    f.write("\n")
    f.write("SofteningTypeOfPartType0   0\n")
    f.write("SofteningTypeOfPartType1   1\n")
    f.write("SofteningTypeOfPartType2   2\n")
    f.write("SofteningTypeOfPartType3   3\n")
    f.write("SofteningTypeOfPartType4   1\n")
    f.write("SofteningTypeOfPartType5   1\n")
    f.write("\n")
    f.write("%MinimumComovingHydroSoftening  1.0\n")
    f.write("%AdaptiveHydroSofteningSpacing  1.2\n")
    f.write("\n")
    f.write("%----- Mesh regularization options\n")
    f.write("CellShapingSpeed       0.5\n")
    f.write("CellShapingFactor      0.8\n")
    f.write("%CellMaxAngleFactor    2.25\n")
    f.write("\n")
    f.write("%----- Subfind\n")
    f.write("ErrTolThetaSubfind     0.7\n")
    f.write("DesLinkNgb             20\n")
    f.write("\n")
    f.write("%----- Black holes\n")
    f.write("%BlackHoleAccretionFactor         100.0\n")
    f.write("%BlackHoleFeedbackFactor          0.07\n")
    f.write("%BlackHoleEddingtonFactor         1.0\n")
    f.write("%SeedBlackHoleMass                1.e-5\n")
    f.write("%MinFoFMassForNewSeed             5.0\n")
    f.write("%DesNumNgbBlackHole               384\n")
    f.write("%BlackHoleMaxAccretionRadius      1.0e5\n")
    f.write("%BlackHoleRadiativeEfficiency     0.2\n")
    f.write("%BHFrictionCoefficient            1.4\n")
    f.write("%BHFrictionAvgTime                0.003\n")
    f.write("\n")
    f.write("%----- Radio Mode\n")
    f.write("%RadioModeMachnumber                0.0075\n")
    f.write("%RadioRelativeBubbleSize            0.1\n")
    f.write("%RadioRelativeBubbleEnergy          0.05\n")
    f.write("%RadioRelativeMaxDist               0.8\n")
    f.write("%RadioModeMetallicityInSolar        1.0\n")
    f.write("\n")
    f.write("%----- Cooling & self-shielding\n")
    f.write("%TreecoolFile                "+arepofilespath+"/TREECOOL_fg_dec11\n")
    f.write("%SelfShieldingFile           "+arepofilespath+"/SelfShielding_Rahmati12\n")
    f.write("%TreecoolFileAGN             "+arepofilespath+"/TREECOOL_AGN\n")
    f.write("%SelfShieldingDensity        0.1295\n")
    f.write("%ObscurationFactor           0.3\n")
    f.write("%ObscurationSlope            0.07\n")
    f.write("%MinMetalTemp                1.0e4\n")
    f.write("%CoolingTablePath            "+arepofilespath+"/Arepo_GFM_Tables/Cooling/cooling_metal_AGN_Compton_self_shielding_Rahmati12.hdf5\n")
    f.write("\n")
    f.write("%----- Enrichment\n")
    f.write("%IMF_MinMass_Msun           0.1\n")
    f.write("%IMF_MaxMass_Msun           100.0\n")
    f.write("%AGB_MassTransferOn         1\n")
    f.write("%SNIa_MassTransferOn        1\n")
    f.write("%SNII_MassTransferOn        1\n")
    f.write("%SNII_MinMass_Msun          8.0\n")
    f.write("%SNII_MaxMass_Msun          100.0\n")
    f.write("%SNIa_Rate_TAU              0.04\n")
    f.write("%SNIa_Rate_Norm             1.3e-3\n")
    f.write("%YieldTablePath             "+arepofilespath+"/Arepo_GFM_Tables/Yields/\n")
    f.write("%DesNumNgbEnrichment          64\n")
    f.write("%MaxNumNgbDeviationEnrichment 1\n")
    f.write("\n")
    f.write("%----- Stellar photometry\n")
    f.write("%PhotometricsTablePath      "+arepofilespath+"/Arepo_GFM_Tables/Photometrics/\n")
    f.write("\n")
    f.write("%----- Refinement\n")
    f.write("%ReferenceGasPartMass    0\n")
    f.write("%TargetGasMassFactor     1.0\n")
    f.write("%RefinementCriterion     1\n")
    f.write("%DerefinementCriterion   1\n")
    f.write("\n")
    f.write("%----- SF & ISM\n")
    f.write("%CritPhysDensity      0       % critical physical density for star formation (in cm^(-3))\n")
    f.write("%MaxSfrTimescale      2.27    % in internal time units (1.5)\n")
    f.write("%CritOverDensity      57.7    % overdensity threshold value\n")
    f.write("%TempSupernova        5.73e7  % in Kelvin (1.0e8)\n")
    f.write("%TempClouds           1000.0  % in Kelvin\n")
    f.write("%FactorEVP            573.0   % (1000.0)\n")
    f.write("%TemperatureThresh    0\n")
    f.write("\n")
    f.write("%----- Winds\n")
    f.write("%WindEnergyIn1e51erg          1.6944001\n")
    f.write("%ThermalWindFactor            3.0\n")
    f.write("%VariableWindVelFactor        3.4641\n")
    f.write("%VariableWindSpecMomentum     0.0\n")
    f.write("%WindFreeTravelMaxTimeFactor  0.025\n")
    f.write("%WindFreeTravelDensFac        0.05\n")
    f.write("%TimeBetOnTheFlyFoF           1.03\n")
    f.write("%MinWindVel                   0.0\n")
    f.write("%WindDumpFactor               0.6\n")
    f.close()
    
def make_arepo_parameter_file_gas(runpath,job_name,run_amd,arepofilespath):
    f = open(runpath + "/param.txt",'w')

    f.write("%----  Relevant files \n")
    f.write("InitCondFile        ./ics\n")
    f.write("OutputDir           ./outputs/\n")
    f.write("SnapshotFileBase    snap\n")
    f.write("OutputListFilename  ExpansionList\n")
    f.write("\n")
    f.write("%---- File formats\n")
    f.write("ICFormat           3 \n")
    f.write("SnapFormat         3 \n")
    f.write("\n")
    f.write("%---- CPU-time limits\n")
    f.write("TimeLimitCPU              2592000   % in seconds\n")
    f.write("CpuTimeBetRestartFile      10800    % in seconds\n")
    f.write("ResubmitOn        0\n")
    f.write("ResubmitCommand   my-scriptfile \n")
    f.write("\n")
    f.write("%----- Memory alloction\n")
    if run_amd:
        f.write("MaxMemSize        3800 \n")
    else:
        f.write("MaxMemSize        1000 \n")
    f.write("\n")
    f.write("%---- Caracteristics of run\n")
    f.write("TimeBegin           0.0078125      % Begin of the simulation\n")
    f.write("TimeMax             1.0            % End of the simulation\n")
    f.write("\n")
    f.write("%---- Basic code options that set the type of simulation\n")
    f.write("ComovingIntegrationOn    1 \n")
    f.write("PeriodicBoundariesOn     1\n")
    f.write("CoolingOn        1\n")
    f.write("StarformationOn      1 \n")
    f.write("\n")
    f.write("%---- Cosmological parameters\n")
    f.write("Omega0                0.3175\n")
    f.write("OmegaLambda           0.6825\n")
    f.write("OmegaBaryon           0.049\n")
    f.write("HubbleParam           0.6711\n")
    f.write("BoxSize               100000.0\n")
    f.write("\n")
    f.write("%---- Output frequency and output paramaters\n")
    f.write("OutputListOn              1 \n")
    f.write("TimeBetSnapshot           0.0\n")
    f.write("TimeOfFirstSnapshot       0.0\n")
    f.write("TimeBetStatistics         0.01\n")
    f.write("NumFilesPerSnapshot       8 \n")
    f.write("NumFilesWrittenInParallel 8\n")
    f.write("\n")
    f.write("%---- Accuracy of time integration\n")
    f.write("TypeOfTimestepCriterion  0 \n")
    f.write("ErrTolIntAccuracy        0.012 \n")
    f.write("CourantFac               0.3 \n")
    f.write("MaxSizeTimestep          0.0025\n")
    f.write("MinSizeTimestep          0.0 \n")
    f.write("\n")
    f.write("%---- Treatment of empty space and temperature limits\n")
    f.write("InitGasTemp            350.0\n")
    f.write("MinGasTemp             15.0\n")
    f.write("MinimumDensityOnStartUp               1.0e-20\n")
    f.write("LimitUBelowThisDensity                0.0\n")
    f.write("LimitUBelowCertainDensityToThisValue  0.0\n")
    f.write("MinEgySpec             0\n")
    f.write("\n")
    f.write("%---- Tree algorithm, force accuracy, domain update frequency\n")
    f.write("TypeOfOpeningCriterion            1\n")
    f.write("ErrTolTheta                       0.7\n")
    f.write("ErrTolForceAcc                    0.0025\n")
    f.write("MultipleDomains                   8 \n")
    f.write("TopNodeFactor                     5\n")
    f.write("ActivePartFracForNewDomainDecomp  0.01\n")
    f.write("\n")
    f.write("%---- Initial density estimate\n")
    f.write("DesNumNgb              64\n")
    f.write("MaxNumNgbDeviation     1 \n")
    f.write("\n")
    f.write("%---- System of units\n")
    f.write("UnitLength_in_cm         3.085678e21        ;  1.0 Kpc \n")
    f.write("UnitMass_in_g            1.989e43           ;  1.0e10 solar masses\n")
    f.write("UnitVelocity_in_cm_per_s 1e5                ;  1 km/sec\n")
    f.write("GravityConstantInternal  0\n")
    f.write("\n")
    f.write("%---- Gravitational softening lengths\n")
    f.write("GasSoftFactor      1.5\n")

    if "LX11" in runpath:
        f.write("SofteningComovingType0    0.610352\n")
        f.write("SofteningComovingType1    0.610352\n")
        f.write("SofteningComovingType2    2.441406\n")
        f.write("SofteningComovingType3    4.882813\n")
        f.write("SofteningComovingType4    195.313\n")
        f.write("SofteningComovingType5    390.625\n")

        f.write("SofteningMaxPhysType0     0.610352\n")
        f.write("SofteningMaxPhysType1     0.610352\n")
        f.write("SofteningMaxPhysType2     2.441406\n")
        f.write("SofteningMaxPhysType3     4.882813\n")
        f.write("SofteningMaxPhysType4     19.5313\n")
        f.write("SofteningMaxPhysType5     39.0625\n")

    if "LX12" in runpath:
        f.write("SofteningComovingType0    0.305175\n")
        f.write("SofteningComovingType1    0.305175\n")
        f.write("SofteningComovingType2    1.220703\n")
        f.write("SofteningComovingType3    2.441406\n")
        f.write("SofteningComovingType4    9.76563\n")
        f.write("SofteningComovingType5    19.5313\n")

        f.write("SofteningMaxPhysType0     0.305175\n")
        f.write("SofteningMaxPhysType1     0.305175\n")
        f.write("SofteningMaxPhysType2     1.220703\n")
        f.write("SofteningMaxPhysType3     2.441406\n")
        f.write("SofteningMaxPhysType4     9.76563\n")
        f.write("SofteningMaxPhysType5     19.5313\n")

    if "LX13" in runpath:
        f.write("SofteningComovingType0    0.15258789\n")
        f.write("SofteningComovingType1    0.15258789\n")
        f.write("SofteningComovingType2    0.61035156\n")
        f.write("SofteningComovingType3    1.22070313\n")
        f.write("SofteningComovingType4    4.88281\n")
        f.write("SofteningComovingType5    9.76563\n")
        
        f.write("SofteningMaxPhysType0     0.15258789\n")
        f.write("SofteningMaxPhysType1     0.15258789\n")
        f.write("SofteningMaxPhysType2     0.61035156\n")
        f.write("SofteningMaxPhysType3     1.22070313\n")
        f.write("SofteningMaxPhysType4     4.88281\n")
        f.write("SofteningMaxPhysType5     9.76563\n")

    if "LX14" in runpath:
        f.write("SofteningComovingType0    0.0762939453125\n")
        f.write("SofteningComovingType1    0.0762939453125\n")
        f.write("SofteningComovingType2    0.30517578\n")
        f.write("SofteningComovingType3    0.61035156\n")
        f.write("SofteningComovingType4    2.44141\n")
        f.write("SofteningComovingType5    4.88281\n")
        
        f.write("SofteningMaxPhysType0     0.0762939453125\n")
        f.write("SofteningMaxPhysType1     0.0762939453125\n")
        f.write("SofteningMaxPhysType2     0.30517578\n")
        f.write("SofteningMaxPhysType3     0.61035156\n")
        f.write("SofteningMaxPhysType4     2.44141\n")
        f.write("SofteningMaxPhysType5     4.88281\n")

    f.write("\n")
    f.write("SofteningTypeOfPartType0   0\n")
    f.write("SofteningTypeOfPartType1   1\n")
    f.write("SofteningTypeOfPartType2   2\n")
    f.write("SofteningTypeOfPartType3   3\n")
    f.write("SofteningTypeOfPartType4   1\n")
    f.write("SofteningTypeOfPartType5   1\n")
    f.write("\n")
    f.write("MinimumComovingHydroSoftening  1.0\n")
    f.write("AdaptiveHydroSofteningSpacing  1.2\n")
    f.write("\n")
    f.write("%----- Mesh regularization options\n")
    f.write("CellShapingSpeed       0.5               \n")
    f.write("CellMaxAngleFactor     2.25\n")
    f.write("\n")
    f.write("%----- Subfind\n")
    f.write("ErrTolThetaSubfind     0.7\n")
    f.write("DesLinkNgb             20\n")
    f.write("\n")
    f.write("%----- Black holes\n")
    f.write("BlackHoleAccretionFactor         100.0\n")
    f.write("BlackHoleFeedbackFactor          0.07\n")
    f.write("BlackHoleEddingtonFactor         1.0\n")
    f.write("SeedBlackHoleMass                1.e-5\n")
    f.write("MinFoFMassForNewSeed             5.0\n")
    f.write("DesNumNgbBlackHole               384\n")
    f.write("BlackHoleMaxAccretionRadius      1.0e5\n")
    f.write("BlackHoleRadiativeEfficiency     0.2\n")
    f.write("BHFrictionCoefficient            1.4\n")
    f.write("BHFrictionAvgTime                0.003\n")
    f.write("\n")
    f.write("%----- Radio Mode\n")
    f.write("RadioModeMachnumber                0.0075\n")
    f.write("RadioRelativeBubbleSize            0.1\n")
    f.write("RadioRelativeBubbleEnergy          0.05\n")
    f.write("RadioRelativeMaxDist               0.8\n")
    f.write("RadioModeMetallicityInSolar        1.0\n")
    f.write("\n")
    f.write("%----- Cooling & self-shielding\n")
    f.write("TreecoolFile                "+arepofilespath+"/TREECOOL_fg_dec11\n")
    f.write("SelfShieldingFile           "+arepofilespath+"/SelfShielding_Rahmati12\n")
    f.write("TreecoolFileAGN             "+arepofilespath+"/TREECOOL_AGN\n")
    f.write("SelfShieldingDensity        0.1295\n")
    f.write("ObscurationFactor           0.3\n")
    f.write("ObscurationSlope            0.07\n")
    f.write("MinMetalTemp                1.0e4\n")
    f.write("CoolingTablePath            "+arepofilespath+"/Arepo_GFM_Tables/Cooling/cooling_metal_AGN_Compton_self_shielding_Rahmati12.hdf5\n")
    f.write("\n")
    f.write("%----- Enrichment\n")
    f.write("IMF_MinMass_Msun           0.1\n")
    f.write("IMF_MaxMass_Msun           100.0\n")
    f.write("AGB_MassTransferOn         1\n")
    f.write("SNIa_MassTransferOn        1\n")
    f.write("SNII_MassTransferOn        1\n")
    f.write("SNII_MinMass_Msun          8.0\n")
    f.write("SNII_MaxMass_Msun          100.0\n")
    f.write("SNIa_Rate_TAU              0.04\n")
    f.write("SNIa_Rate_Norm             1.3e-3\n")
    f.write("YieldTablePath             "+arepofilespath+"/Arepo_GFM_Tables/Yields/\n")
    f.write("DesNumNgbEnrichment          64\n")
    f.write("MaxNumNgbDeviationEnrichment 1\n")
    f.write("\n")
    f.write("%----- Stellar photometry\n")
    f.write("PhotometricsTablePath      "+arepofilespath+"/Arepo_GFM_Tables/Photometrics/\n")
    f.write("\n")
    f.write("%----- Refinement\n")
    f.write("ReferenceGasPartMass    0\n")
    f.write("TargetGasMassFactor     1.0\n")
    f.write("RefinementCriterion     1\n")
    f.write("DerefinementCriterion   1\n")
    f.write("\n")
    f.write("%----- SF & ISM\n")
    f.write("CritPhysDensity      0       % critical physical density for star formation (in cm^(-3))\n")
    f.write("MaxSfrTimescale      2.27    % in internal time units (1.5)\n")
    f.write("CritOverDensity      57.7    % overdensity threshold value\n")
    f.write("TempSupernova        5.73e7  % in Kelvin (1.0e8)\n")
    f.write("TempClouds           1000.0  % in Kelvin\n")
    f.write("FactorEVP            573.0   % (1000.0)\n")
    f.write("TemperatureThresh    0\n")
    f.write("\n")
    f.write("%----- Winds\n")
    f.write("WindEnergyIn1e51erg          1.6944001\n")
    f.write("ThermalWindFactor            3.0\n")
    f.write("VariableWindVelFactor        3.4641\n")
    f.write("VariableWindSpecMomentum     0.0\n")
    f.write("WindFreeTravelMaxTimeFactor  0.025\n")
    f.write("WindFreeTravelDensFac        0.05\n")
    f.write("TimeBetOnTheFlyFoF           1.03\n")
    f.write("MinWindVel                   0.0\n")
    f.write("WindDumpFactor               0.6\n")
    f.close()


def run_arepo(suite_paths,arepo_file_path,lx_list,submit=True,dmonly=True):
    job_name_list = []
    current_jobs,jobids,jobstatus,namd,nhyper = getcurrentjobs()
    tallyhyper = nhyper
    tallyamd = namd
    
    for folder in suite_paths:
        topfolder = "/".join(folder.split("/")[:8])
        #print topfolder
        folder_single = folder.split("/")[-1]
        halo_label = folder_single.split("_")[0]
        nrvir_label = folder_single.split("NV")[1][0]
        lx_label = folder_single.split("LX")[1][1:2]

        if dmonly:
            job_name = "d"+halo_label+"X"+lx_label+"N"+nrvir_label+folder_single.split(halo_label+"_")[1][:2]
        else:
            job_name = "g"+halo_label+"X"+lx_label+"N"+nrvir_label+folder_single.split(halo_label+"_")[1][:2]

        lastsnap = '098'
      
        if os.path.isfile(folder + "/ics.0.hdf5") and os.path.getsize(folder + "/ics.0.hdf5") > 0 and \
            job_name not in current_jobs and job_name not in job_name_list:
            if not os.path.isdir(folder + "/outputs/snapdir_"+lastsnap) and lx_label in lx_list:
                print
                print "COPYING AREPO FILES...",folder
                mkdir_outputs = "mkdir -p " + folder + "/outputs/"
                
                if "LX11" in folder or "LX12" in folder:
                    bakpath = folder+"/outputs/restartfiles/restart.*"
                    nbak = len(glob.glob(bakpath))
                    #print nbak
                    if nbak == 64:
                        ncores = 64
                        pmgrid = 512
                        queue = "AMD64"
                        run_amd = True
                        tallyamd += 1
                    elif nbak == 120:
                        queue = "HyperNodes"
                        ncores = 120
                        pmgrid = 512
                        run_amd = False
                        tallyhyper += 1
                    elif nbak == 0:
                        if tallyamd <= tallyhyper:
                            ncores = 64
                            pmgrid = 512
                            queue = "AMD64"
                            run_amd = True
                            tallyamd += 1
                        if tallyamd > tallyhyper:
                            queue = "HyperNodes"
                            ncores = 120
                            pmgrid = 512
                            run_amd = False
                            tallyhyper += 1
                    
                    #print ncores
                    print "RUNNING AMD:",run_amd
                    print "%i (%i) JOBS ON HYPER (AMD64) NODES!" % (nhyper,namd)
                    #sys.exit()

                if "LX13" in folder:
                    ncores = 64
                    pmgrid = 1024
                    queue = "AMD64"
                    run_amd = True
                    tallyamd += 1

                if "contamination" in folder:
                    file_times  = "cp /bigbang/data/bgriffen/forpaul/exe/ExpansionList_contam " + folder + "/ExpansionList"
                else:
                   file_times  = "cp /bigbang/data/bgriffen/forpaul/exe/ExpansionList_merge " + folder + "/ExpansionList"

                if len(glob.glob(folder+"/outputs/snapdir*")) > 0 and len(glob.glob(folder+"/outputs/restartfiles/*.bak")) == ncores:
                    restart_flag = True
                else:
                    restart_flag = False

                #print folder

                if len(glob.glob(folder+"/outputs/snapdir*")) == 0 and nbak == 0:
                    print "OUTPUTS/ IS CLEAN > STARTING AT START!"
                    skip_flag = False
                elif len(glob.glob(folder+"/outputs/snapdir*")) == 0 and nbak == ncores:
                    print "OUTPUTS/ IS CLEAN > STARTING AT START!"
                    skip_flag = False
                elif len(glob.glob(folder+"/outputs/snapdir*")) == 0 and nbak != ncores:
                    print "OUTPUTS/ IS CLEAN > STARTING AT START!"
                    skip_flag = True
                elif len(glob.glob(folder+"/outputs/snapdir*")) != 0 and nbak == ncores:
                    print "RESTARTING FROM RESTART FILES > " + glob.glob(folder+"/outputs/snapdir*")[-1].split("/")[-1]
                    skip_flag = False
                elif len(glob.glob(folder+"/outputs/snapdir*")) != 0 and nbak != ncores:
                    print "WANTED TO RESTART FROM RESTART FILES BUT # RESTART FILES != ncores > " + glob.glob(folder+"/outputs/snapdir*")[-1].split("/")[-1]
                    skip_flag = True
                    print "# restart files:",len(glob.glob(folder+"/outputs/restartfiles/*.bak")),"# cores:",ncores

                if not skip_flag:
                    #file_param  = "cp " + gadget_file_path + "/arepo/param_1"+lx_label+"_" + queue + ".txt " + folder + "/param.txt"
                    if dmonly:
                        file_exe    = "cp /bigbang/data/bgriffen/forpaul/exe/Arepo_"+str(pmgrid)+ "_dm "  + folder + "/Arepo"
                        file_config = "cp /bigbang/data/bgriffen/forpaul/exe/Config_"+str(pmgrid)+"_dm.sh " + folder + "/Config.sh"
                        make_arepo_parameter_file_dm(folder,job_name,run_amd,arepofilespath="/bigbang/data/bgriffen/lib/arepofiles/")
                    else:
                        file_exe    = "cp /bigbang/data/bgriffen/forpaul/exe/Arepo_"+str(pmgrid)+ "_gas "  + folder + "/Arepo"
                        file_config = "cp /bigbang/data/bgriffen/forpaul/exe/Config_"+str(pmgrid)+"_gas.sh " + folder + "/Config.sh"
                        make_arepo_parameter_file_gas(folder,job_name,run_amd,arepofilespath="/bigbang/data/bgriffen/lib/arepofiles/")
                    
                    cmd_copy_all_files = [mkdir_outputs,file_times,file_exe,file_config]
                    subprocess.call([";".join(cmd_copy_all_files)],shell=True)
                    #make_gadget_submission_script(folder,job_name,restart_flag,gadget3,run_amd)
                    make_arepo_submission_script(folder,job_name,restart_flag,run_amd)
                    
                    cd_folder = "cd " + folder
                    cmd_submit_arepo = "sbatch sarepo"

                    if submit == True:
                        print "SUBMITTING AREPO..."
                        subprocess.call([cd_folder+"; "+cmd_submit_arepo],shell=True)
                        job_name_list.append(job_name)

                else:
                    print "NOT SUBMITTING!"

def getcentext(filename):
    with open(filename, 'r') as fp:
        lines = []
        for i in xrange(6):
            lines.append(fp.readline().strip('#\n'))

    #print lines
    lines = map(float, lines)
    return tuple(lines)
    
def run_music(suite_paths,music_path,lagr_path,lx_list,dmonly=True):
     for folder in suite_paths:
         
         folder_single = folder.split("/")[-1]
         halo_label = folder_single.split("_")[0]
         nrvir_label = folder_single.split("NV")[1][0]
         lx_label = folder_single.split("LX")[1][1:2]
         #job_name = halo_label[:4]+"X"+lx_label+"N"+nrvir_label+folder_single.split("_")[-1]
         if dmonly: 
            job_name = "dI"+halo_label[1:]+"X"+lx_label+"N"+nrvir_label+folder_single.split(halo_label+"_")[1][:2]
         else:
            job_name = "gI"+halo_label[1:]+"X"+lx_label+"N"+nrvir_label+folder_single.split(halo_label+"_")[1][:2]

         job_name_list = []
         current_jobs,jobids,jobstatus,namd,nhyper = getcurrentjobs()

         skip = False
         runsingle = False
         if os.path.isfile(folder + "/ics.0.hdf5"):
            if os.path.getsize(folder + "/ics.0.hdf5") > 0:
                skip = True
                print "skipping",folder

         #print skip
         if not skip and lx_label in lx_list:
             #print folder
             if job_name not in current_jobs and job_name not in job_name_list:
                 print
                 print "RUNNING:",folder_single
                 print "MAKING MUSIC SUBMISSION SCRIPT..."
    
                 make_music_submission_script(folder,folder_single,job_name,runsingle=runsingle)
                 
                 print "COPYING MUSIC FILES..."
                 if runsingle:
                     print "running in SINGLE PRECISION!"
                     print folder
                     cmd_cp_music = "cp " + music_path + "_singleprec " + folder
                     cmd_rm_music = "rm " + folder + "/MUSIC"
                     subprocess.call([cmd_rm_music],shell=True)
                 else:
                     cmd_cp_music = "cp " + music_path + " " + folder

                 subprocess.call([cmd_cp_music],shell=True)
             
                 print "CONSTRUCTING MUSIC CONFIGURATION FILES..."
                 lagr_file = lagr_path + folder_single.split("_")[0]
                 master_music_cfg_dest = folder + "/" + folder_single + ".conf"
                 
                 region_point_file = lagr_path+"H"+halo_label[1:]+"NRVIR"+nrvir_label
                 print region_point_file
                 seed = int(halo_label[1:])
                 
                 if os.path.isfile(region_point_file + ".head"):
                     if halo_label+"_B" in folder_single:
                        ictype = 'box'
                        
                     if halo_label+"_E" in folder_single:
                        ictype = "ellipsoid"
    
                     if halo_label+"_C" in folder_single:
                        ictype = "convex_hull"
    
                 #print master_music_cfg_dest,folder_single,ictype,seed,region_point_file,nrvir_label
                 #sys.exit()
                 #print ictype
                 make_music_file(master_music_cfg_dest,ictype,seed,region_point_file,nrvir_label,runsingle)
    
                 #np.loadtxt(lagr_path+"H"+halo_label[1:]+".head")
                 #with open(lagr_path+"H"+halo_label[1:]+".head") as myfile:
    
                 #make_LX11_musicfile(master_music_cfg_dest,ictype,seed,region_point_file)

                 print "SUBMITTING INITIAL CONDITIONS..."
                 cd_folder = "cd " + folder
                 cmd_submit_ics = "sbatch smusic"
                 subprocess.call([cd_folder+"; "+cmd_submit_ics],shell=True)
                 job_name_list.append(job_name)
                 #sys.exit()

def run_music_higher_levels(halo_geometries,base_path,music_path,lagr_path,run_halos,lx_list,dmonly):
    for hid in run_halos:
        halo_name = "H"+str(hid)
        for LX in ["11","12","13"]:
            if LX[1] in lx_list: # and halo_name not in skip_list:
                #hpath = htils.hid_hpath_lx(hid,int(LX))
                #print hpath
                geometry = halo_geometries["H"+str(hid)].replace("_","")[:2]

                #geometry,lxdump,nvir = htils.get_zoom_params(hpath)
                nvir = str(3)
                lxuse = str(LX)
                #if dmonly: 
                #    folder = base_path + halo_name + "/d" + halo_name +  "_"+geometry+"_Z127_P7_LN7_LX"+lxuse+"_O4_NV"+nvir 
                #else:
                folder = base_path + halo_name + "/" + halo_name +  "_"+geometry+"_Z127_P7_LN7_LX"+lxuse+"_O4_NV"+nvir 
                if not os.path.isdir(folder):
                    #print halo_name
                    cmd_make_next_level = "mkdir -p " + folder
                    print cmd_make_next_level
                    subprocess.call([cmd_make_next_level],shell=True)
                    
                sub_suite_paths = glob.glob(base_path + "/"+halo_name+"/H*")
                run_music(sub_suite_paths,music_path,lagr_path,lx_list,dmonly=dmonly)
    
def make_music_file(master_music_cfg_dest,boxtype,seed,region_point_file,nrvir_label,runsingle):
    haloid = float(master_music_cfg_dest.split("/")[-1].split("_")[0][1:])
    folder_single = master_music_cfg_dest.split("/")[-1].replace("cfg","")
    halo_label = folder_single.split("_")[0]
    lx_label = folder_single.split("LX")[1][1:2]

    f = open(master_music_cfg_dest,'w')
    f.write("[setup]\n")
    f.write("boxlength            = 100\n")
    f.write("zstart               = 127\n")
    f.write("levelmin             = 7\n")
    f.write("levelmin_TF          = 10\n")
    f.write("levelmax             = 1"+lx_label+"\n")
    f.write("padding              = 7\n")
    f.write("overlap              = 4\n")
    f.write("region               = " + boxtype + "\n")
    
    if boxtype == 'box':
        centx,centy,centz,dx,dy,dz = getcentext(region_point_file+".head")
        
        if "_BA_" in master_music_cfg_dest:
            extx=dx
            exty=dy
            extz=dz
        if "_BB_" in master_music_cfg_dest:
            extx=dx*1.2 
            exty=dy*1.2 
            extz=dz*1.2 
        if "_BC_" in master_music_cfg_dest:
            extx=dx*1.4 
            exty=dy*1.4 
            extz=dz*1.4 
        if "_BD_" in master_music_cfg_dest:
            extx=dx*1.6 
            exty=dy*1.6 
            extz=dz*1.6 

        f.write("ref_center           = " + str(centx) + "," + str(centy) + "," + str(centz)+"\n")
        f.write("ref_extent           = " + str(extx) + "," + str(exty) + "," + str(extz)+"\n")
    
    if boxtype == 'ellipsoid':
        #print nrvir_label
        if "_EA_" in master_music_cfg_dest:
            hipadding = 1.0
        if "_EB_" in master_music_cfg_dest:
            hipadding = 1.1
        if "_EC_" in master_music_cfg_dest:
            hipadding = 1.2      
        if "_ED_" in master_music_cfg_dest:
            hipadding = 1.3  
        if "_EX_" in master_music_cfg_dest:
            hipadding = 1.05    

        f.write("hipadding            = " + str(hipadding) + "\n")

    f.write("region_point_file    = " + region_point_file + "\n")
    f.write("align_top            = no\n")

    f.write("baryons              = no\n")

    f.write("use_2LPT             = no\n")
    f.write("use_2LLA             = no\n")
    f.write("periodic_TF          = yes\n")
    f.write("\n")
    f.write("[cosmology]\n")
    f.write("Omega_m              = 0.3175\n")
    f.write("Omega_L              = 0.6825\n")
    f.write("Omega_b              = 0.049\n")

    f.write("H0                   = 67.11\n")
    f.write("sigma_8              = 0.8344\n")
    f.write("nspec                = 0.9624\n")
    f.write("transfer             = eisenstein\n")
    f.write("\n")
    f.write("[random]\n")
    f.write("seed[10]              = 34567\n")

    for LX in xrange(1,int(lx_label)+1):
        seed = int(haloid*(float(LX % 10)))
        f.write("seed[1"+str(LX)+"]              = " + str(seed) + "\n")
        #print LX,seed

    f.write("\n")
    f.write("[output]\n")
    f.write("format               = arepo\n")
    f.write("filename             = ics.hdf5\n")
    f.write("arepo_num_files      = 8\n")
    f.write("arepo_doubleprec     = 1\n")
    f.write("\n")
    f.write("[poisson]\n")
    f.write("fft_fine             = yes\n")
    f.write("accuracy             = 1e-05\n")
    f.write("pre_smooth           = 3\n")
    f.write("post_smooth          = 3\n")
    f.write("smoother             = gs\n")
    f.write("laplace_order        = 6\n")
    f.write("grad_order           = 6\n")
    f.close()

def getcurrentjobs():
    pipemyq = 'squeue -o "%i,%j,%t,%P" -u bgriffen > currentqueue.out'
    subprocess.call(';'.join([pipemyq]),shell=True)
    lines = [line.strip() for line in open('currentqueue.out')]
    subprocess.call(["rm currentqueue.out"],shell=True)
    currentjobs = []
    jobstatus = []
    jobids = []
    queuename = []
    #print lines
    for i in xrange(1,len(lines)):
        linesplit = lines[i].split(",")
        jobids.append(linesplit[0])
        currentjobs.append(linesplit[1])
        jobstatus.append(linesplit[2])
        queuename.append(linesplit[3])

    #print queuename
    namds = queuename.count('AMD64')
    nhyper = queuename.count('HyperNodes')

    return currentjobs,jobids,jobstatus,namds,nhyper

def make_top_halo_folders(lagr_path):
    for lagrfile in glob.glob(lagr_path+"H*.head"):
        cmd_make_folder = "mkdir -p /bigbang/data/bgriffen/forpaul/halos/" + lagrfile.split("/")[-1].replace(".head","").split("NRVIR")[0]
        subprocess.call([cmd_make_folder],shell=True)

def make_destination_folders_clean(base_path,suite_names,lx,nrvir):
    for folder in glob.glob(base_path+"H*"):
        haloid = folder.split("halos/")[-1].split("_")[0]
        new_folder_name = haloid + "_BB_Z127_P7_LN7_LX"+str(lx)+"_O4_NV"+str(nrvir)
        #folder_single = folder_path.split("/")[-1]
        for suite in suite_names:
            cmd_make_folder = "mkdir -p " + folder +  "/contamination_suite/" + new_folder_name + "_" + suite 
            subprocess.call([cmd_make_folder],shell=True)

def make_destination_folders(base_path,suite_names,lx,nrvir):
    for folder in glob.glob(base_path+"H*"):
        haloid = folder.split("halos/")[-1].split("_")[0]
        #folder_single = folder.split("/")[-1]
        for suite in suite_names:
            new_folder_name = haloid + "_" + suite + "_Z127_P7_LN7_LX"+str(lx)+"_O4_NV"+str(nrvir)
            #old_folder_name = haloid + "_BB_Z127_P7_LN7_LX"+str(lx)+"_O4_NV"+str(nrvir)
            #print new_folder_name
            #cmd_rm_folder = "rm -rf " + folder +  "/contamination_suite/" + old_folder_name + "_" + suite 
            #subprocess.call([cmd_rm_folder],shell=True)
            cmd_make_folder = "mkdir -p " + folder + "/contamination_suite/" + new_folder_name
            subprocess.call([cmd_make_folder],shell=True)

