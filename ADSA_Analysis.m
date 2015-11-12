classdef ADSA_Analysis < handle

% Analysis class for a 3-dimensional framed structure
    
    % Private properties go here
    properties (Access = private)
            nnodes
            Coordinates
            ConcenLoads
            Fixity
            nele
            Ends
            A
            Izz
            Iyy
            J
            Zzz
            Zyy
            Ayy
            Azz
            E
            v
            webdir
            DistribLoads
            Nodes
            Elements
            StructureStiffnessMatrix
            FeF
            P
    end
    
    % Public methods go here
    methods (Access = public)
        %% Constructor
 
        function self = ADSA_Analysis(nnodes, coord, concen, fixity, nele, ends,A,Izz,Iyy,J,Zzz,Zyy,...
		Ayy,Azz,E,v,webdir, w)
            self.nnodes=nnodes;
            self.Coordinates=coord;
            self.ConcenLoads=concen;
            self.Fixity=fixity;
            self.nele=nele;
            self.Ends=ends;
            self.A=A;
            self.Izz=Izz;
            self.Iyy=Iyy;
            self.J=J;
            self.Zzz=Zzz;
            self.Zyy=Zyy;
            self.Ayy=Ayy;
            self.Azz=Azz;
            self.E=E;
            self.v=v;
            self.webdir=webdir;
            self.DistribLoads=w;
            %Creating a vector of Node Objects
            self.Nodes=CreateNodes(self);
            %Creating a vector of Element Objects
            self.Elements=CreateElements(self);
            %Assembling the structure stiffness matrix from the global
            %stiffness matrix of each element
            self.StructureStiffnessMatrix = CreateStiffnessMatrix(self);
            %Assembling the fixed end forces and the concentrated loads on
            %the structure from those acting on each element
            [self.FeF, self.P] = CreateLoadVectors(self);
        end
        
        %% Run the analysis
        function [AFLAG, DEFL, REACT, ELE_FOR]=RunAnalysis(self)
        	[AFLAG, DEFL, REACT, ELE_FOR]=GetMastan2Returns(self)
            [Error]=ComputeError(self)
        end
        
        function [AFLAG, DEFL, REACT, ELE_FOR]=GetMastan2Returns(self)
            %Calls ComputeDisplacementReactions
            %RecoverElementForces
            ELE_FOR=RecoverElementForces(self, DEFL);
            %CheckKffMatrix
        end
        
        
        
    end
    
    % Private methods go here
    methods (Access = private)
        
              
        function error=ComputeError(self)
            %Space for computing
        end
        
        function [DEFL, REACT]=ComputeDisplacementReactions(self)
            
            %Classify DOF
            [freeDOF, fixedDOF, knownDOF]= ClassifyDOF(self);
            
            %CreatStiffnessMatrix
            [Kff, Kfn, Knf, Knn, Ksf, Ksn]= ComputeStiffnessSubMatrices(self, freeDOF, fixedDOF, knownDOF);
           
            %ComputStiffnessSubMatrices
            %CreateLoadVector
        end
        
        %Method to create a vector of node objects
        function Nodes=CreateNodes(self)
            %Loop through all nodes
            for i=1:self.nnodes
                %Create Node objects and store in vector called Nodes
                Nodes(i)=ADSA_Node(self.Coordinates(i,:), i);
            end
        end
            
        %Method to create a vector of element objects using the node
        %objects and element properties
        function Elements=CreateElements(self)
            %Loop through all elements
            for i=1:self.nele
                %Store each object in a vector with index corresponding to
                %its number.  Send the node objects that correspond to each
                %end of the element and the element properties as arguments
                Elements(i)=ADSA_Element([self.Nodes(self.Ends(i, 1)), ...
                    self.Nodes(self.Ends(i,2))],self.A(i),...
                    self.Izz(i),self.Iyy(i),self.J(i),self.Cw(i),...
                    self.Zzz(i),self.Zyy(i),self.Ayy(i),self.Azz(i),...
                    self.E(i),self.v(i),self.Fy(i),self.YldSurf(i,:),...
                    self.Wt(i),self.webdir(i,:), self.DistribLoads(i,:));
            end
        end
        
        
        function StructureStiffnessMatrix = CreateStiffnessMatrix(self)
            
            %Initializing the K matrix to all zeros
            K = zeros(self.nnodes*6); 
            
            % Looping over all the element objects
            for i= 1:self.nele
                %Assembling the structure stiffness matrix by adding the
                %values of the global stiffness matrix of the elements to
                %the corresponding DOFs in the structure stiffness matrix
                K(GetElementDOF(self.Elements(i)),GetElementDOF(self.Elements(i))) = ...
                    K(GetElementDOF(self.Elements(i)),GetElementDOF(self.Elements(i)))... 
                    + GetGlobalStiffness(self.Elements(i));
            end
            disp('Structure Stiffness Matrix= ');
            disp(K);
            K = sparse(K);
            StructureStiffnessMatrix = K;
        end
        
        function [Kff, Kfn, Knf, Knn, Ksf, Ksn]= ComputeStiffnessSubMatrices(self, freeDOF, fixedDOF, knownDOF)
                      
            Kff= K(freeDOF, freeDOF);
            Kfn= K(freeDOF, knownDOF);
            Knf= K(knownDOF, freeDOF);
            Knn= K(knownDOF, knownDOF);
            Ksf= K(fixedDOF, freeDOF);
            Ksn= K(fixedDOF, knownDOF);
        end
            
            
        function [FeF, P] = CreateLoadVectors(self)
            
            %Initializing the FeF matrix to all zeros
            FeF = zeros(self.nnodes*6,1);
            
            % Looping over all the element objects
            for i= 1:self.nele
                %Assembling the matrix of the fixed end forces of the 
                %structure by adding the values of the global FeF matrix of
                %the elements to the corresponding DOFs in the structure 
                %FeF matrix
                FeF(GetElementDOF(self.Elements(i)))=FeF(GetElementDOF(self.Elements(i)))...
                    + GetFixedEndForcesGlobal(self.Elements(i));
            end
            
            %Initializing the P matrix to all zeros
            P = zeros(self.nnodes*6,1);
            
            % Looping over all the element objects
            for i= 1:self.nele
                %Assembling the matrix of the concentrated forces on the 
                %structure by adding the values of the concentrated forces
                %on a specific DOF to the corresponding DOF in the P matrix
                %of the structure
                P(GetElementDOF(self.Elements(i))) = P(GetElementDOF(self.Elements(i)))...
                 + [self.ConcenLoads(GetNodeDOF(self.Nodes(self.Ends(i, 1))));...
                 self.ConcenLoads(GetNodeDOF(self.Nodes(self.Ends(i, 2))))];
            end
            disp('FeF= ');
            disp(FeF);
            disp('Concentrated Loads= ');
            disp(P);
        end
        
    function [freeDOF, fixedDOF, knownDOF]= ClassifyDOF(self)
        
        %Transpose the fixity property matrix so linear indexing matches
        %the number of the DOF of each node
        fixitytrans=self.fixity';
        
        %Find indices of free DOF's
        freeDOF=find(isnan(fixitytrans));
        
        %Find indices of fixed (supported) DOF's
        fixedDOF=find(fixitytrans==0);
        
        %Find indices of DOF's with assigned displacement values
        knownDOF=find(fixitytrans~=0 & ~isnan(fixitytrans));
        
    end
    
    function AFLAG=CheckKffMatrix(self)
        
        %Obtain Kff matrix
        Kff=ComputeStiffnessSubMatrices(self)
        
        %For efficiency, generate an estimate of the condition number
        kappa=condest(self.Kff);
        fprintf('Condition number is %.5f', kappa)
        
        %Estimate the number of significant digits will be lost
        %Estimate that P-S=log(kappa) where P is input significant digits
        %and S is number of reliable significant digits returned
        lostDigits=log(kappa);
        fprintf('Results likely to loose %d significant digits', lostDigits)
        
    end
    
    
    function [DEFL, REACT]=ComputeDisplacementsReactions(self)
        
        %Space for computing
    end
    
    function ELE_FOR=RecoverElementForces(self, DEFL)
        
        %Initialize element forces to a matrix of zeroes
        ELE_FOR=zeros(self.nele, 12);
        
        %Loop over all elements
            for i=1:self.nele
                ELE_FOR(i, 1:12)=ComputeForces(self.Elements(i), ...
                    [DEFL(self.ends(i,1), 1:6),DEFL(self.ends(i,2), 1:6)]);
            end
    end
      
    end
    
end 
    