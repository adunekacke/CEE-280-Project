classdef ADSA_Analysis < handle

% Analysis class for a 3-dimensional framed structure
    
    % Private properties go here
    properties (Access = private)
            nnodes
            Coordinates
            concen
            fixity
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
           
    end
    
    % Public methods go here
    methods (Access = public)
        %% Constructor
 
        function self = ADSA_Analysis(nnodes, coord, concen, fixity, nele, ends,A,Izz,Iyy,J,Zzz,Zyy,...
		Ayy,Azz,E,v,webdir, w)
            self.nnodes=nnodes;
            self.Coordinates=coord;
            self.concen=concen;
            self.fixity=fixity;
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

        end
        
        %% Run the analysis
        function [AFLAG, DEFL, REACT, ELE_FOR]=RunAnalysis(self)
            %Call function which calculates values Mastan is expecting
        	[AFLAG, DEFL, REACT, ELE_FOR]=GetMastan2Returns(self);
            
            %Call function to compute the percent error in loads at free
            %DOF's
            Error=ComputeError(self, DEFL);
            
            %Display to Command Window

            disp('Error in Loads at Free DOF''s:')
            fprintf('%.16f\n\n', Error);

        end
        
        function [AFLAG, DEFL, REACT, ELE_FOR]=GetMastan2Returns(self)
            %Calls ComputeDisplacementReactions
            [DEFL, REACT]=ComputeDisplacementReactions(self);
            
            %RecoverElementForces
            ELE_FOR=RecoverElementForces(self, DEFL);
            
            %CheckKffMatrix
            AFLAG= CheckKffMatrix(self);
        end
        
        
        
    end
    
    % Private methods go here
    methods (Access = private)
        
              
        function error=ComputeError(self, DEFL)
                        
            %Get DOF's that are free and known
            [freeDOF, fixedDOF, knownDOF]=ClassifyDOF(self);
            
            %Transpose DEFL for use of linear indexing
            DEFL=DEFL';

            %Get vectors of free and known deflections
            deltaf=DEFL(freeDOF);
            deltan=DEFL(knownDOF);

            %Get stiffness submatrices
            [Kff, Kfn]= ComputeStiffnessSubMatrices(self, freeDOF, fixedDOF, knownDOF);
            
            %Calculate Loads at free DOF's
            Pfback=Kff*deltaf+Kfn*deltan;
            
            %Create actual load vector
            [Pfreal, ~, ~, FEFreal] = CreateLoadVectors(self, freeDOF, fixedDOF, knownDOF);
            realLoads=Pfreal-FEFreal;

            %Compute error in loads at free DOF's
            error=realLoads-Pfback;
                       
        end
        
        function [DEFL, REACT]=ComputeDisplacementReactions(self)
            
            %Classify DOF
            [freeDOF, fixedDOF, knownDOF]= ClassifyDOF(self);

            %Obtaining the Stiffness Sub Matrices
            [Kff, Kfn, Knf, Knn, Ksf, Ksn]= ComputeStiffnessSubMatrices(self, freeDOF, fixedDOF, knownDOF);
            
            %Obtaining the applied loads, displacements and fixed end
            %forces
            [Pf, Ps, Pn, FeFf, FeFs, FeFn, DeltaN] = CreateLoadVectors(self, freeDOF, fixedDOF, knownDOF);
            
            % Displacements at the Free DOFs
            DeltaF= Kff\(Pf - FeFf - Kfn*DeltaN);
            
            %Forces/Reactions at Fixed DOFs
            Ps= Ksf*DeltaF + Ksn*DeltaN + FeFs - Ps;
            
            %Forces/Reactions at Known DOFs
            Pn= Knf*DeltaF + Knn*DeltaN + FeFn - Pn;
            
            %Initializing DEFL to zeros of 6 rows (No. of DOFs in a node)
            %and total number of nodes as number of columns
            DEFL= zeros(6, self.nnodes);
            
            %Placing the values of DeltaF in the indeces corresponding to
            %the Free DOFs
            DEFL(freeDOF)= DeltaF;
            
            %Placing the values of DeltaN in the indeces corresponding to
            %the Known DOFs
            DEFL(knownDOF)= DeltaN;
            
            %Placing Zero as values in the indeces corresponding to
            %the Fixed DOFs as there will be zero displacements there
            DEFL(fixedDOF)= 0;
            
            %Transposing the DEFL matrix to have 6 columns (corresponding 
            %to 6 DOFs in a node) and number of nodes as the number of
            %rows- The format which MASTAN2 requires
            DEFL= DEFL';
            
            %Initializing REACT to zeros of 6 rows (No. of DOFs in a node)
            %and total number of nodes as number of columns
            REACT= zeros(6, self.nnodes);
            
            %Placing Zero as values in the indeces corresponding to
            %the Free DOFs as there will be No reactions there
            REACT(freeDOF)= 0;
            
            %Placing the values of Pn in the indeces corresponding to
            %the Known DOFs
            REACT(knownDOF)= Pn;
            
            %Placing the values of Ps in the indeces corresponding to
            %the Fixed DOFs
            REACT(fixedDOF)= Ps;
            
            %Transposing the REACT matrix to have 6 columns (corresponding 
            %to 6 DOFs in a node) and number of nodes as the number of
            %rows- The format which MASTAN2 requires
            REACT= REACT';
                

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
                    self.Izz(i),self.Iyy(i),self.J(i),...
                    self.Zzz(i),self.Zyy(i),self.Ayy(i),self.Azz(i),...
                    self.E(i),self.v(i), self.webdir(i,:), ...
                    self.DistribLoads(i,:));
            end
        end
        
        %Method to generate the complete structural stiffness matrix by
        %combining global stiffness matrices of individual elements
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
           
            %Store as a sparse matrix
            K = sparse(K);
            StructureStiffnessMatrix = K;
        end
        
        %Method to split the structural stiffness matrix into the funtional
        %pieces needed.
        function [Kff, Kfn, Knf, Knn, Ksf, Ksn]= ComputeStiffnessSubMatrices(self, freeDOF, fixedDOF, knownDOF)

            K= self.StructureStiffnessMatrix;

            Kff= K(freeDOF, freeDOF);
            Kfn= K(freeDOF, knownDOF);
            Knf= K(knownDOF, freeDOF);
            Knn= K(knownDOF, knownDOF);
            Ksf= K(fixedDOF, freeDOF);
            Ksn= K(fixedDOF, knownDOF);
        end
        
            
        %Method to create the load vectors from the concentrated loads and 
        %distributed loads and divide into pieces that correspond to 
        %degrees of freedom with free, fixed, and known displacements.
        function [Pf, Ps, Pn, FeFf, FeFs, FeFn, DeltaN] = CreateLoadVectors(self, freeDOF, fixedDOF, knownDOF)
            
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
            
            concen_t= self.concen';
            
            Pf= concen_t(freeDOF);
            Ps= concen_t(fixedDOF);
            Pn= concen_t(knownDOF);
            
            FeFf= FeF(freeDOF);
            FeFs= FeF(fixedDOF);
            FeFn= FeF(knownDOF);
            
            
            fixitytrans= self.fixity';
            DeltaN= fixitytrans(knownDOF);

        end
        
    function [freeDOF, fixedDOF, knownDOF]= ClassifyDOF(self)
        
        %Transpose the fixity property matrix so linear indexing matches
        %the number of the DOF of each node
        fixitytrans= self.fixity';
        
        %Find indices of free DOF's
        freeDOF= find(isnan(fixitytrans));
        
        %Find indices of fixed (supported) DOF's
        fixedDOF= find(fixitytrans==0);
        
        %Find indices of DOF's with assigned displacement values
        knownDOF= find(fixitytrans~=0 & ~isnan(fixitytrans));
        
    end
    
    %Method to check the conditioning of K 
    function AFLAG= CheckKffMatrix(self)
        
        %Get DOF's
        [freeDOF, fixedDOF, knownDOF]= ClassifyDOF(self);
        
        %Obtain Kff matrix
        Kff=ComputeStiffnessSubMatrices(self, freeDOF, fixedDOF, knownDOF);
        
        %For efficiency, generate an estimate of the condition number
        kappa= condest(Kff);
        fprintf('Condition number is %d\n\n', kappa)
        
        %Estimate the number of significant digits will be lost
        %Estimate that P-S=log(kappa) where P is input significant digits
        %and S is number of reliable significant digits returned
        lostDigits=log10(kappa);
        fprintf('Results likely to lose %.1f significant digits\n\n', lostDigits)
        
        %Assign AFLAG based on how many significant digits will be lost
        if lostDigits > 13
            AFLAG= 0;
        else
            AFLAG= 1;
        end
    end
    
    
    function [DEFL, REACT]=ComputeDisplacementsReactions(self)
        
        %Space for computing
    end
    
    %Method to recover element forces from the element method and store in
    %the Mastan format
    function ELE_FOR=RecoverElementForces(self, DEFL)
        
        %Initialize element forces to a matrix of zeroes
        ELE_FOR=zeros(self.nele, 12);
        
        %Loop over all elements, calling public method of element and
        %passing element object and deflections 
            for i=1:self.nele
                ELE_FOR(i, 1:12)=ComputeForces(self.Elements(i), ...
                    [DEFL(self.Ends(i,1), 1:6),DEFL(self.Ends(i,2), 1:6)]);
            end
    end
    
    end
    
end 
    
