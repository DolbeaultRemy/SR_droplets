classdef PhysicsConstants < handle
    properties (Constant)
        % CODATA
        PlanckConstant=6.62607015E-34;
        PlanckConstantReduced=6.62607015E-34/(2*pi);
        FineStructureConstant=7.2973525698E-3;
        ElectronMass=9.10938291E-31;
        GravitationalConstant=6.67384E-11;
        ProtonMass=1.672621777E-27;
        AtomicMassUnit=1.660539066E-27; 
        BohrRadius=5.2917721067E-11; 
        BohrMagneton=9.274009994E-24; 
        BoltzmannConstant=1.38064852E-23;
        StandardGravityAcceleration=9.80665;
        SpeedOfLight=299792458;
        StefanBoltzmannConstant=5.670373E-8;
        ElectronCharge=1.602176634E-19;
        VacuumPermeability=1.25663706212E-6; 
        DielectricConstant=8.8541878128E-12;
        ElectronGyromagneticFactor=-2.00231930436153;
        AvogadroConstant=6.02214076E23;
        ZeroKelvin = 273.15;
        GravitationalAcceleration = 9.80553;
        
        % Dy specific constants
        Dy164Mass                    = 163.929174751*1.660539066E-27;
        Dy164IsotopicAbundance       = 0.2826;
        DyMagneticMoment             = 9.93*9.274009994E-24;
    end
    
    methods
        function pc = PhysicsConstants()
        end
    end
    
end
