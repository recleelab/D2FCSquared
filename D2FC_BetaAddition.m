function [model,Data] = ApoptosisSimBiologyModelFunction(Scenario,LocOfData,TimeIntervalSeconds)

%==========================================================================
%                        Model Initiation 
%==========================================================================
model = sbiomodel("NFkB Model");


%==========================================================================
%                        Model Compartments
%==========================================================================
Cytoplasm = addcompartment(model,'Cytoplasm',"Capacity",3.33);
Nucleus   = addcompartment(model,'Nucleus',"Capacity",1.0);


helper = Helper(model,"micromole","second");




%==========================================================================
%                Apply Sum of Gaussians for IKK Trajectory 
%==========================================================================
helper.Scalar("a1_coeff",0);
helper.Scalar("b1_coeff",0);
helper.Scalar("c1_coeff",0);
helper.Scalar("a2_coeff",0);
helper.Scalar("b2_coeff",0);
helper.Scalar("c2_coeff",0);
helper.Scalar("a3_coeff",0);
helper.Scalar("b3_coeff",0);
helper.Scalar("c3_coeff",0);
helper.Scalar("a4_coeff",0);
helper.Scalar("b4_coeff",0);
helper.Scalar("c4_coeff",0);


%Sum of Four Guassians for the IKK 
helper.Scalar("IKKSpots",0,false)

sumOfFourGaussians = ['a1_coeff*exp(-((time-b1_coeff)/c1_coeff)^2)' ...
                     '+a2_coeff*exp(-((time-b2_coeff)/c2_coeff)^2)' ...
                     '+a3_coeff*exp(-((time-b3_coeff)/c3_coeff)^2)' ...
                     '+a4_coeff*exp(-((time-b4_coeff)/c4_coeff)^2)'];

addrule(model, "IKKSpots = "+sumOfFourGaussians, 'repeatedAssignment');




%==========================================================================
%                           Model Species
%==========================================================================
helper.AddSpecies(Cytoplasm,"NFkB",0)
helper.AddSpecies(Cytoplasm,"IkBa",0)
helper.AddSpecies(Cytoplasm,"IkBaNFkB",0.3)

helper.AddSpecies(Nucleus,"NFkB",0)
helper.AddSpecies(Nucleus,"IkBa",0)
helper.AddSpecies(Nucleus,"IkBaNFkB",0)

helper.AddSpecies(Cytoplasm,"pIkBa",0)
helper.AddSpecies(Cytoplasm,"pIkBaNFkB",0)


helper.AddSpecies(Cytoplasm,"IKK",0)
helper.AddSpecies(Cytoplasm,"IKKi",0)
helper.AddSpecies(Cytoplasm,"IKKn",0.8)

helper.AddSpecies(Cytoplasm,"tIkBa",0)
helper.AddSpecies(Cytoplasm,"tCompetitor",0)
helper.AddSpecies(Cytoplasm,"Competitor",0)

helper.AddSpecies(Cytoplasm,"tA20",0)
helper.AddSpecies(Cytoplasm,"A20",0)

helper.AddSpecies(Cytoplasm,"IkBb_t",0)
helper.AddSpecies(Cytoplasm,"IkBb",0)
helper.AddSpecies(Nucleus,"IkBb",0)
helper.AddSpecies(Nucleus,"IkBbNFkB",0)
helper.AddSpecies(Cytoplasm,"IkBbNFkB",0)
helper.AddSpecies(Cytoplasm,"pIkBb",0)
helper.AddSpecies(Cytoplasm,"pIkBbNFkB",0)

%==========================================================================
%                            Parameters
%==========================================================================

parameters_d2fc;


%%%%%%%%%%%% d2fc options %%%%%%%%%%%%%%%%

parameters_d2fc;


helper.Scalar("TR",0,false)

%Kv is the ratio of cytoplasm to nuclear area and needs to stay 
%constant. 
helper.Scalar("kv",kv);


%==========================================================================
%                        Model Reactions
%==========================================================================
%Binding and unbinding of NFkB and IkB in Cytoplasm
helper.SecondOrder('ka1a',ka1a);
helper.FirstOrder("kd1a",kd1a);
helper.AddMassActionEquation( ...
           'Cytoplasm.IkBa + Cytoplasm.NFkB <-> Cytoplasm.IkBaNFkB', ...
            {'ka1a','kd1a'})

%Reaction rates of NFkB transferring between cytoplasm and nucleus while
%accounting for change of volume 
helper.FirstOrder("ki1",ki1);
helper.FirstOrder("ke1",ke1);

[r2,r3,r4,r5]= GenerateCytoplasmNucleurTransfer(model,"NFkB","ki1","ke1","kv");

%Binding and unbinding of NFkB and IkB in Nucleus
helper.Scalar("f_ka1a",0.1)
helper.Scalar("f_kd1a",0.1)

helper.SecondOrder("ka1a_nucleus",1E-14)
helper.FirstOrder("kd1a_nucleus",1E-14)

addrule(model,"ka1a_nucleus = f_ka1a* ka1a","initialAssignment");
addrule(model,"kd1a_nucleus = f_kd1a* kd1a","initialAssignment");
helper.AddMassActionEquation( ...
            'Nucleus.IkBa + Nucleus.NFkB <-> Nucleus.IkBaNFkB', ...
             {'ka1a_nucleus','kd1a_nucleus'})

%Decay of IkBa in the nucleus and cytplasm
helper.FirstOrder("c4a",c4a)
helper.AddMassActionEquation("Nucleus.IkBa -> null",{'c4a'})
helper.AddMassActionEquation("Cytoplasm.IkBa -> null",{'c4a'})

%Exit of NFkB/IkB complex from the nucelus to the Cytoplasm 
helper.FirstOrder("ke2a",ke2a);
r7 = RateExitNucleus(model,"IkBaNFkB","ke2a","kv");

%Reaction rates for IkB transferring between the cytoplasm and nucleus
%while accounting for change of volume 
helper.FirstOrder("ki3a",ki3a);
helper.FirstOrder("ke3a",ke3a);
[r8,r9,r10,r11] = GenerateCytoplasmNucleurTransfer(model,"IkBa","ki3a","ke3a","kv");

%IkB Cytoplasmic Degragation on NFkB. NFkB is generated 
helper.FirstOrder("c5a",c5a);
helper.AddMassActionEquation( ...
                             "Cytoplasm.IkBaNFkB -> Cytoplasm.NFkB", ...
                             {'c5a'})

%Neutral IKK (IKKn) converting to active IKK (IKK) 
helper.FirstOrder("ka",1E-3);
helper.AddReaction( ...
         "Cytoplasm.IKKn -> Cytoplasm.IKK",...
         "TR*IKKSpots*ka*Cytoplasm.IKKn")

%IKK converting to inactive IKKn
helper.FirstOrder("ki",ki)
helper.AddMassActionEquation("Cytoplasm.IKK -> Cytoplasm.IKKi",{'ki'})

%Inactive IKK (IKKn) converting to neutural IKK (IKKn)
helper.FirstOrder("kp",kp)
helper.Concentration("kbA20",kbA20)
helper.AddReaction("Cytoplasm.IKKi -> Cytoplasm.IKKn", ...
    "kp*Cytoplasm.IKKi*kbA20/(kbA20+A20)")


%Formation of IkBa mRNA transcript in response to NFkB in the nucleus 
helper.Scalar("c1a",0.5);
helper.Scalar("h",h);
helper.Concentration("k",k);
helper.ZeroOrder("rs_a",3.0800e-06)
RateOfIkBaTranscription = "rs_a.*(1 + c1a.*(Nucleus.NFkB./k).^h ./( (Nucleus.NFkB./k).^h + 1))";
helper.AddReaction("null -> Cytoplasm.tIkBa",RateOfIkBaTranscription)

%Decay of IkBa mRNA transcript. Follows basic mass action kinetics 
helper.FirstOrder("c3a",c3a);
helper.AddMassActionEquation("Cytoplasm.tIkBa -> null",{'c3a'})

%Transcript of IkBa can make IkBa in the cytoplasm 
helper.FirstOrder("c2a",c2a)
helper.AddReaction("null -> Cytoplasm.IkBa",...
                   "c2a*Cytoplasm.tIkBa")


%Competitor Transcript Formation 
helper.Concentration("k4",k4)
helper.AddReaction("null -> Cytoplasm.tCompetitor",...
 "rs_a* c1a* (Nucleus.NFkB/k)^(h+1) /( (Nucleus.NFkB/k)^(h+1) + (Competitor/k4)^(h+1)  + 1)")

%Competitor Transcript Decay 
helper.FirstOrder("c6a",c6a)
helper.AddMassActionEquation("Cytoplasm.tCompetitor -> null",{'c6a'})

%Competitor Protein Creation 
helper.AddReaction("null -> Cytoplasm.Competitor", ... 
              "c2a*Cytoplasm.tCompetitor")

%Competitor Protein Decay 
helper.AddMassActionEquation("Cytoplasm.Competitor -> null",{'c4a'})

%A20 transcript formation [tA20]
helper.ZeroOrder("c1",c1)
helper.Concentration("k2",k2)
helper.AddReaction("null -> Cytoplasm.tA20",...
 "c1* (Nucleus.NFkB/k)^(h+1) /( (Nucleus.NFkB/k)^(h+1) + (Cytoplasm.Competitor/k2)^(h+1)  + 1)")

%A20 Trasncript Decay 
helper.FirstOrder("c3",c3)
helper.AddMassActionEquation("Cytoplasm.tA20 -> null",{'c3'})

%A20 Creation 
helper.FirstOrder("c2",c2)
helper.AddReaction("null -> Cytoplasm.A20", ... 
              "c2*Cytoplasm.tA20")

%A20 Protein Decay 
helper.FirstOrder("c4",c4)
helper.AddMassActionEquation("Cytoplasm.A20 -> null",{'c4'})

%Phosphorlating IkBa 
helper.SecondOrder("kc1a",kc1a)
helper.AddReaction("Cytoplasm.IkBa -> Cytoplasm.pIkBa", ...
                  "kc1a*Cytoplasm.IKK*Cytoplasm.IkBa")

%Decay of pIkBa 
helper.FirstOrder("kt1a",kt1a)
helper.AddMassActionEquation("Cytoplasm.pIkBa -> null",{'kt1a'})

%Phosphorlating pIkBaNFkB 
helper.SecondOrder("kc2a",kc2a)
helper.AddReaction("Cytoplasm.IkBaNFkB -> Cytoplasm.pIkBaNFkB", ...
                  "kc2a*Cytoplasm.IKK*Cytoplasm.IkBaNFkB")

%Phosphorlated pIkBaNFkB will break apart where IkBa will decay and NFkB
%will be released 
helper.FirstOrder("kt2a",kt2a)
helper.AddMassActionEquation("Cytoplasm.pIkBaNFkB -> Cytoplasm.NFkB", ...
                            {'kt2a'})

%==========================================================================
%IkBb addition was based off of Hoffmann et al. 2002 (The IkB-NF-kB
%Signaling Module: Temporal Control and Selective Activation) and Werner et
%al. 2005 (Stimulus Specificity of Gene Expression Programs Determined by
%Temporal Control of IKK activity). Parameters from Werner et al used. 
%==========================================================================
%IkBb Transcript Formation 
helper.ZeroOrder("rs_b",7.12e-7)
helper.AddMassActionEquation("null -> Cytoplasm.IkBb_t",{'rs_b'})

helper.FirstOrder("rd_b",2.8e-4) %Same between Werner and Hoffman 
helper.AddMassActionEquation("Cytoplasm.IkBb_t -> null",{'rd_b'})

%IKBb Translation and degragation 
helper.FirstOrder("ps_c_b",4.08e-3) %Same between Werner and Hoffman 
helper.AddReaction("null -> Cytoplasm.IkBb","ps_c_b*Cytoplasm.IkBb_t")

helper.FirstOrder("pd_n_b",0.003)
helper.AddMassActionEquation("Cytoplasm.IkBb -> null",{'pd_n_b'})

%Transport of IkBb in and out of the nucleus
helper.FirstOrder("in_b",3e-4)
helper.FirstOrder("ex_b",2e-4)
GenerateCytoplasmNucleurTransfer(model,"IkBb","in_b","ex_b","kv");

%Binding and unbinding of NFkB and IkBb in Cytoplasm and nucleus
helper.SecondOrder("a_c_bn",0.5)
helper.FirstOrder("d_c_bn",kd1a)
helper.AddMassActionEquation( ...
            'Cytoplasm.IkBb + Cytoplasm.NFkB <-> Cytoplasm.IkBbNFkB', ...
             {'a_c_bn','d_c_bn'})

helper.SecondOrder("a_n_bn",0.5)
helper.FirstOrder("d_n_bn",1e-6)

helper.Scalar("f_a_n_bn",0.1)
helper.Scalar("f_d_n_bn",0.1)

addrule(model,"a_n_bn = f_a_n_bn* a_c_bn","initialAssignment");
addrule(model,"d_n_bn = f_d_n_bn* d_c_bn","initialAssignment");

helper.AddMassActionEquation( ...
            'Nucleus.IkBb + Nucleus.NFkB <-> Nucleus.IkBbNFkB', ...
             {'a_n_bn','d_n_bn'})

%Degragation of IkBb bound to NFkB complex 
helper.FirstOrder("pd_c_2bn",1e-6)
helper.AddMassActionEquation("Cytoplasm.IkBbNFkB -> Cytoplasm.NFkB",{'pd_c_2bn'})

%Export of IkBbNFkB complexes in the nucleus to the cytoplasm 
helper.FirstOrder("ex_2bn",7.1e-3)
RateExitNucleus(model,"IkBbNFkB","ex_2bn","kv");


%Phosphorlating IkBa 
helper.SecondOrder("kc1b",kc1a)
helper.AddReaction("Cytoplasm.IkBb -> Cytoplasm.pIkBb", ...
                  "kc1b*Cytoplasm.IKK*Cytoplasm.IkBb")

%Decay of pIkBa 
helper.FirstOrder("kt1b",kt1a)
helper.AddMassActionEquation("Cytoplasm.pIkBb -> null",{'kt1b'})

%Phosphorlating pIkBaNFkB 
helper.SecondOrder("kc2b",kc2a)
helper.AddReaction("Cytoplasm.IkBbNFkB -> Cytoplasm.pIkBbNFkB", ...
                  "kc2b*Cytoplasm.IKK*Cytoplasm.IkBbNFkB")

%Phosphorlated pIkBaNFkB will break apart where IkBa will decay and NFkB
%will be released 
helper.FirstOrder("kt2b",kt2a)
helper.AddMassActionEquation("Cytoplasm.pIkBbNFkB -> Cytoplasm.NFkB", ...
                            {'kt2b'})
%==========================================================================
%Observables for model output 
%==========================================================================
addobservable(model,"TotalNFkB_Cyt","Cytoplasm.NFkB + Cytoplasm.IkBaNFkB + Cytoplasm.pIkBaNFkB+Cytoplasm.pIkBbNFkB+Cytoplasm.IkBbNFkB");
addobservable(model,"TotalNFkB_Nuc","Nucleus.NFkB +Nucleus.IkBaNFkB+  Nucleus.IkBbNFkB");
addobservable(model,"FractionNuc2Cyt","(TotalNFkB_Nuc./kv)./TotalNFkB_Cyt");
addobservable(model,"RelativeNFkB_Nuc","TotalNFkB_Nuc./TotalNFkB_Nuc(1)");

addobservable(model,"Total_IkBalphaCytoplasm","Cytoplasm.IkBa + Cytoplasm.IkBaNFkB");
addobservable(model,"Total_IkBalphaPhoshorlated","Cytoplasm.pIkBa + Cytoplasm.pIkBaNFkB");


addobservable(model,"Total_IkBalphaNucleus","Nucleus.IkBa + Nucleus.IkBaNFkB");
addobservable(model,"FractionIkBaCytoplasm","Total_IkBalphaCytoplasm./(Total_IkBalphaCytoplasm) ");
addobservable(model,"FractionIkBaNucleus","Total_IkBalphaNucleus./(Total_IkBalphaNucleus) ");



%Need to Log Parameter States 
configsetObj = getconfigset(model);
% Get the current StatesToLog
currentStatesToLog = configsetObj.RuntimeOptions.StatesToLog;

pars2Log = {'PioneeringScale','KDNA', 'nDNA','c1a','k','h','kv','NFkBStability','nUnbinding','kd_bind','nCooperativity'};
for i = 1:length(pars2Log)
    ithParameter = pars2Log{i};
    ithParameter = sbioselect(model,"Type","parameter","Name",ithParameter);
    currentStatesToLog = [currentStatesToLog; ithParameter];
end 
 configsetObj.RuntimeOptions.StatesToLog = currentStatesToLog;

end 

