classdef ADSA_Analysis < handle

%Analysis class for a three-dimentional elastic frame structure
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Functions Called
%     Public:
%       ADSA_Analysis      Constructs the analysis object and stores its
%                            properties
%       RunAnalysis        Conducts analysis of the structure, performing
%                            stifffness method computations and returning 
%                            the structure's deflections, reactions, and
%                            element forces to Mastan2 
%       GetMastan2Returns  Getter function to get deflections, reactions,
%                            and element forces again after analysis has 
%                            been run inititially
%     Private:
%       CreateNodes        Loops through all the nodes and creates the 
%                            node objects
%       CreateElements     Loops through all the elements and creates the 
%                            element objects
%       ComputeStiffnessSubMatrices   Gets the stiffness matrices for each
%                            element, sums coinciding terms, and extracts
%                            the necessary sub-matrices
%       ClassifyDOF        Classifies the degrees of freedom as free,
%                            fixed, or known based on the fixity matrix
%       ComputeError       Back-calculates the applied forces from the 
%                            deflections and compares to actual loads
%       CheckKff           Checks the condition number and estimates lost
%                            digits to determine if analysis will be 
%                            successful
%       ComputeDisplacementReactions    Computes the displacements and 
%                            reactions using the stiffness method
%       CreateLoadVectors  Separates concentrated loads and fixed end
%                            forces based on degree of freedom 
%                            classification once it gets the fixed end 
%                            forces from the element class
%       RecoverElementForces  Using the displacements, call the element 
%                            method to compute internal forces 
%       
%
%  Dictionary of Variables
%     Input Information and Properties:
%       nnodes         ==  total number of nodes
%       coordinates(i,1:3)   ==  node i's coordinates
%                            coord(i,1) = X coordinate
%                            coord(i,2) = Y coordinate
%                            coord(i,3) = Z coordinate
%       concen(i,1:6)  ==  concentrated loads for node i's 6 d.o.f.
%                            concen(i,1) = force in global X direction
%                            concen(i,2) = force in global Y direction
%                            concen(i,3) = force in global Z direction
%                            concen(i,4) = moment about global X axis
%                            concen(i,5) = moment about global Y axis
%                            concen(i,6) = moment about global Z axis
%       fixity(i,1:6)  ==  prescribed displacements for node i's 6 d.o.f.
%                          Note: A free d.o.f. will have a value of NaN
%                          and hence, you will find the Matlab function
%                          isnan very useful.
%                          Examples: If fixity(15,3) is set to NaN, then node 15's
%                                      Z-disp component is free;
%                                    If fixity(2,6) is set to 0.0, then node 2's
%                                      Z-rotation component is supported;
%                                    If fixity(5,2) is set to -2.1, then node 5's
%                                      Y-disp component is supported and defined
%                                      with a settlement of -2.1 units.
%                            fixity(i,1) = prescribed disp. in global X direction
%                            fixity(i,2) = prescribed disp. in global Y direction
%                            fixity(i,3) = prescribed disp. in global Z direction
%                            fixity(i,4) = prescribed rotation about global X axis
%                            fixity(i,5) = prescribed rotation about global Y axis
%                            fixity(i,6) = prescribed rotation about global Z axis
%       nele           ==  total number of elements
%       ends(i,1:14)   ==  element i's nodal information
%                            ends(i,1) = start node #
%                            ends(i,2) = finish node #
%                            ends(i,3) = flag to indicate whether or not flexural
%                            moments are released at start node.  ends(i,3)=0 both not
%                            released (rigid connection); ends(i,3)=1 both flexural
%                            moments are released (pinned connection); ends(i,3)=2
%                            at least one of the flexural moments are partially or fully
%                            released (see below for connection stiffness attributes)
%                            ends(i,4) = flag to indicate whether or not flexural
%                            moments are released at finish node.  ends(i,4)=0 both not
%                            released (rigid connection); ends(i,4)=1 both flexural
%                            moments are released (pinned connection); ends(i,4)=2
%                            at least one of the flexural moments are partially or fully
%                            released (see below for connection stiffness attributes)
%                            ends(i,5) = flag to indicate the degree of warping
%                            restraint at start node.  ends(i,5)=0 warping free;
%                            ends(i,5)=1 warping fixed; ends(i,5)=2 warping continuous
%                            ends(i,6) = flag to indicate the degree of warping
%                            restraint at finish node.  ends(i,6)=0 warping free;
%                            ends(i,6)=1 warping fixed; ends(i,6)=2 warping continuous
%                            ends(i,7) = rotational spring stiffness at the start
%                            node and about element i's local z-z axis.
%                            ends(i,8) = rotational spring stiffness at the start
%                            node and about element i's local y-y axis.
%                            ends(i,9) = rotational spring stiffness at the finish
%                            node and about element i's local z-z axis.
%                            ends(i,10) = rotational spring stiffness at the finish
%                            node and about element i's local y-y axis.
%                            ends(i,11) = connection moment capacity Mpz at the start
%                            node and about element i's local z-z axis.
%                            ends(i,12) = connection moment capacity Mpy at the start
%                            node and about element i's local y-y axis.
%                            ends(i,13) = connection moment capacity Mpz at the finish
%                            node and about element i's local z-z axis.
%                            ends(i,14) = connection moment capacity Mpy at the finish
%                            node and about element i's local y-y axis.
%       A(i)           ==  element i's cross sectional area
%       Izz(i)         ==  element i's moment of inertia about its local z-z axis
%       Iyy(i)         ==  element i's moment of inertia about its local y-y axis
%       J(i)           ==  element i's torsional constant
%       Zzz(i)         ==  element i's plastic section modulus about its local z-z axis
%       Zyy(i)         ==  element i's plastic section modulus about its local y-y axis
%       Ayy(i)         ==  element i's effective shear area along its local y-y axis
%       Azz(i)         ==  element i's effective shear area along its local z-z axis
%       E(i)           ==  element i's material elastic modulus, Young's Modulus
%       v(i)           ==  element i's material Poisson's ratio
%       webDir(i,1:3)  ==  element i's unit web vector.  This is a unit vector
%                          that defines the element's local y-y axis with respect
%                          to the global coordinate system.  It is based on the
%                          structure's undeformed geometry.
%                              webdir(i,1) = x component of element's unit web vector
%                              webdir(i,2) = y component of element's unit web vector
%                              webdir(i,3) = z component of element's unit web vector
%                          NOTE: An element's 3x3 rotation matrix, [g], is constructed
%                          as follows: First, calculate a unit vector, x_vect, that
%                          describes the element's local x-axis. Second, take the
%                          cross product of x_vect and webdir(i,:) to obtain z_vect,
%                          i.e. z_vect = cross(x_vect,webdir(i,:)). Third, set z_vect 
%                          to a unit vector, i.e. z_vect = z_vect/norm(z_vect).
%                          Finally, the first row of [g] is x_vect, its second row is
%                          webdir(i,:), and its third row is z_vect.
%       distribLoads(i,1:3)==  element i's uniform load which references its
%                            local coordinate system
%                              w(i,1) = x component of uniform load
%                              w(i,2) = y component of uniform load
%                              w(i,3) = z component of uniform load
%       Nodes             Vector of Node objects of length nnodes
%       Elements          Vector of Element objects of length nele
%       Kff               Matrix of free-free stiffness terms
%       Kfn               Matrix of free-known stiffness terms
%       Knf               Matrix of known-free stiffness terms
%       Knn               Matrix of known-known stiffness terms
%       Ksf               Matrix of fixed-free stiffness terms
%       Ksn               Matrix of fixed-known stiffness terms
%
%     Local Variables:
%       m                Stores row/column length of the Structure
%                        Stiffness Matrix
%       el_DOF           Stores the vector of DOFs of an element while
%                        computing the sparse Structure Stiffness Matrix
%       el_Stiff         Stores the global stiffness matrix of an element
%                        while computing the sparse Structure Stiffness 
%                        Matrix
%       row              Vector consisting of the DOFs corresponding to the
%                        row indeces of the non-zero stiffness values while
%                        computing the sparse Structure Stiffness Matrix
%       col              Vector consisting of the DOFs corresponding to the
%                        column indeces of the non-zero stiffness values 
%                        while computing the sparse Structure Stiffness 
%                        Matrix
%       stiff            Vector Consisting of non-zero stiffness values
%                        corresponding to the DOFs 'row' and 'col'
%       count            An incrementing variable to append values to the
%                        'row', 'col' and 'stiff' vectors
%       K                Complete structural stiffness matrix
%       fixityTrans      The fixity matrix transposed for linear indexing
%       freeDOF          The numbers of the free degrees of freedom
%       fixedDOF         The numbers of the fixed degrees of freedom
%       knownDOF         The numbers of the known degrees of freedom
%       deltaf           The displacements at the free degrees of freedom
%       deltan           Displacements at the known degrees of freedom
%       backPf           The back-calculated loads on the system
%       realPf           The real concentrated loads from the input
%       realFEF          Real fixed end forces from input distributed loads
%       realLoads        Real loads on system (concentrated and distributed)
%       Pf               Concentrated loads at free degrees of freedom
%       Ps               Concentrated loads at fixed supports
%       Pn               Concentrated loads at supports with known disp.
%       FEFf             Fixed end forces at free degrees of freedom
%       FEFs             Fixed end forces at fixed supports
%       FEFn             Fixed end forces at supports w/ known displacment
%       Rs               Reactions at fixed supports
%       Rn               Reactions at supports with known displacements
%       FeF              All fixed end forces obtained from elements
%       concen_t         Transposed matrix of concentrated loads
%
%     Output Information:
%       DEFL(i,1:6)      ==  node i's calculated 6 d.o.f. deflections
%                              DEFL(i,1) = displacement in X direction
%                              DEFL(i,2) = displacement in Y direction
%                              DEFL(i,3) = displacement in Z direction
%                              DEFL(i,4) = rotation about X direction
%                              DEFL(i,5) = rotation about Y direction
%                              DEFL(i,6) = rotation about Z direction
%       REACT(i,1:6)     ==  reactions for supported node i's 6 d.o.f.
%                              REACT(i,1) = force in X direction
%                              REACT(i,2) = force in Y direction
%                              REACT(i,3) = force in Z direction
%                              REACT(i,4) = moment about X direction
%                              REACT(i,5) = moment about Y direction
%                              REACT(i,6) = moment about Z direction
%       ELE_FOR(i,1:12)  ==  element i's internal forces and moments
%                            Note: All values reference the element's local
%                                  coordinate system.
%                              ELE_FOR(i,1)  = x-force at start node
%                              ELE_FOR(i,2)  = y-force at start node
%                              ELE_FOR(i,3)  = z-force at start node
%                              ELE_FOR(i,4)  = x-moment at start node
%                              ELE_FOR(i,5)  = y-moment at start node
%                              ELE_FOR(i,6)  = z-moment at start node
%                              ELE_FOR(i,7)  = x-force at end node
%                              ELE_FOR(i,8)  = y-force at end node
%                              ELE_FOR(i,9)  = z-force at end node
%                              ELE_FOR(i,10) = x-moment at end node
%                              ELE_FOR(i,11) = y-moment at end node
%                              ELE_FOR(i,12) = z-moment at end node
%       AFLAG            ==  logical flag to indicate if a successful
%                            analysis has been completed
%                              AFLAG = 1     Successful
%                              AFLAG = 0     Unstable Structure
%                              AFLAG = inf   No analysis code available
%       kappa            condition number of the free-free stiffness matrix
%       lostDigits       estimate of how many digits will be lost between
%                              the input and the returns
%       error            estimate of the error in the calculations;
%                              difference between input loads and back-
%                              calculated values
%
    
  % Private properties 
    properties (Access = private)
            nnodes
            coordinates
            concen
            fixity
            nele
            ends
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
            webDir
            distribLoads
            Nodes
            Elements
            Kff
            Kfn
            Knf
            Knn
            Ksf
            Ksn
            AFLAG
            DEFL
            REACT
            ELE_FOR
            error
    end
    
    
  % Public methods go here
    methods (Access = public)
        
        %% Constructor
        function self = ADSA_Analysis(nnodes, coordinates, concen, fixity,...
                     nele, ends,A,Izz,Iyy,J,Zzz,Zyy,Ayy,Azz,E,v,webDir,...
                     distribLoads)
            self.nnodes=nnodes;
            self.coordinates=coordinates;
            self.concen=concen;
            self.fixity=fixity;
            self.nele=nele;
            self.ends=ends;
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
            self.webDir=webDir;
            self.distribLoads=distribLoads;
            self.Nodes=CreateNodes(self);
            self.Elements=CreateElements(self);
            [self.Kff, self.Kfn, self.Knf, self.Knn, ...
                 self.Ksf, self.Ksn]= ComputeStiffnessSubMatrices(self);
            
        end
        
        
        %% Run the analysis
        function [AFLAG, DEFL, REACT, ELE_FOR]=RunAnalysis(self)
            
            %CheckKffMatrix
            AFLAG= CheckKffMatrix(self);
            
            %If analyis is succesful, run computations
            if AFLAG==1
                %Calls ComputeDisplacementReactions to get deflections and
                %reactions
                [DEFL, REACT]=ComputeDisplacementReactions(self);

                %RecoverElementForces to compute the internal forces in the
                %elements
                ELE_FOR=RecoverElementForces(self);
                
                %Compute the error in the associated with the loads
                error= ComputeError(self);
                
                %Display error to Command Window
                disp('Error in Loads at Free DOF''s:')
                fprintf('%.16f\n', error);
                
            %If unsuccessful, halt analysis
            else
                self.DEFL= zeros(self.nnodes, 6);
                self.REACT= zeros(self.nnodes, 6);
                self.ELE_FOR= zeros(self.nele, 12);
                disp('Error is not applicable; structure is unstable.')
            end
           
        end
        
        %Method to make the returns that Mastan2 expects
        function [AFLAG, DEFL, REACT, ELE_FOR, error]= GetMastan2Returns(self)
            
            AFLAG= self.AFLAG;
            DEFL= self.DEFL;
            REACT= self.REACT;
            ELE_FOR= self.ELE_FOR;
            error= self.error;
        end
    end
    
  % Private methods go here
    methods (Access = private)
        
        %Method to create a vector of node objects
        function Nodes=CreateNodes(self)
            %Loop through all nodes
            for i=1:self.nnodes
                %Create Node objects and store in vector called Nodes
                Nodes(i)=ADSA_Node(self.coordinates(i,:), i);
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
                Elements(i)=ADSA_Element([self.Nodes(self.ends(i, 1)), ...
                    self.Nodes(self.ends(i,2))],self.A(i),...
                    self.Izz(i),self.Iyy(i),self.J(i),...
                    self.Zzz(i),self.Zyy(i),self.Ayy(i),self.Azz(i),...
                    self.E(i),self.v(i), self.webDir(i,:), ...
                    self.distribLoads(i,:));
            end
        end
        
        
        %Method to Create the structure stiffness matrix and then split it
        %into the funtional pieces needed.
        function [Kff, Kfn, Knf, Knn, Ksf, Ksn]= ...
                                        ComputeStiffnessSubMatrices(self)
            %Creating the structure stiffness matrix and Computing
            %Stiffness Sub Matrices methods have been merged into one.
            %Keeping them as separate methods required the K matrix (in its
            %sparse form) to be stored as a property and the
            %ComputeStiffnessSubMatrices function being called thrice in
            %the analysis. Now only the required 6 out of the 9 sub parts
            %of the K matrix are being stored and the function is also
            %required to be called just once in the analysis.
            
            %Creating a sparse Structure Stiffness Matrix by
            %combining the non-zero values of the global stiffness matrices
            %of individual elements
            
            %Storing the size (rows or columns) of the complete Structure 
            %Stiffness Matrix
            m= self.nnodes*6;
            
            %Initializing 'count' to 1
            count= 1;
            
            % Looping over all the element objects
            for z= 1:self.nele
                
                %Storing the zth Element's DOFs as a vector in 'el_DOF'
                el_DOF= GetElementDOF(self.Elements(z));
                
                %Storing the zth Element's Global Stiffness Matrix as
                %'el_Stiff'
                el_Stiff= GetGlobalStiffness(self.Elements(z));
                
                for x= 1:12
                    for y= 1:12
                        %Concentrating only on non-zero values
                        if el_Stiff(x,y)~=0 
                            
                            % The DOF el_DOF(x) corresponding to non-zero 
                            % stiffness value gets appended to the 'row' 
                            % vector
                            row(count)= el_DOF(x);
                            
                            % The DOF el_DOF(y) corresponding to non-zero 
                            % stiffness value gets appended to the 'col' 
                            % vector
                            col(count)= el_DOF(y);
                            
                            % The value of stiffness corresponding to the
                            % DOFs of el_DOF(x) and el_DOF(y) gets appended
                            % to the 'stiff' vector
                            stiff(count)= el_Stiff(x,y);
                            count= count + 1;
                        end
                    end
                end
            end
            
            % Using the format of the sparse matrix: K= sparse(i,j,v,m,n)
            % the Structure Stiffness Matrix is stored as a sparse matrix
            % Since Total No. of rows= Total No. of columns, m=n= nnodes*6
            K= sparse(row, col, stiff, m, m); 
            
            [freeDOF, fixedDOF, knownDOF]=  ClassifyDOF(self);

            Kff= K(freeDOF, freeDOF);
            Kfn= K(freeDOF, knownDOF);
            Knf= K(knownDOF, freeDOF);
            Knn= K(knownDOF, knownDOF);
            Ksf= K(fixedDOF, freeDOF);
            Ksn= K(fixedDOF, knownDOF);
        end
        
        
        %Method to classify the degrees of freedom: free, fixed, and known
        function [freeDOF, fixedDOF, knownDOF]= ClassifyDOF(self)

            %Transpose the fixity property matrix so linear indexing
            %matches the number of the DOF of each node
            fixityTrans= self.fixity';
        
            %Find indices of free DOF's
            freeDOF= find(isnan(fixityTrans));

            %Find indices of fixed (supported) DOF's
            fixedDOF= find(fixityTrans==0);

            %Find indices of DOF's with assigned displacement values
            knownDOF= find(fixityTrans~=0 & ~isnan(fixityTrans));
        
        end

        
        %Method to compute the error in the analysis due to lost digits
        function error= ComputeError(self)
                        
            %Get DOF's that are free and known
            [freeDOF, fixedDOF, knownDOF]=ClassifyDOF(self);
            
            %Transpose DEFL for use of linear indexing
            DEFL=self.DEFL';

            %Get vectors of free and known deflections
            deltaF=DEFL(freeDOF);
            deltaN=DEFL(knownDOF);

            %Calculate Loads at free DOF's
            backPf=self.Kff*deltaF+self.Kfn*deltaN;
            
            %Create actual load vector
            [realPf, ~, ~, realFEF] = CreateLoadVectors(self, freeDOF, ...
                                                      fixedDOF, knownDOF);
            realLoads= realPf-realFEF;

            %Compute error in loads at free DOF's
            self.error= realLoads-backPf;
            error= self.error;
                       
        end
        
        
        %Method to check the conditioning of K 
        function AFLAG= CheckKffMatrix(self)
        
            %For efficiency, generate an estimate of the condition number
            kappa= condest(self.Kff);
            fprintf('Condition number is %d\n\n', kappa)
        
            %Estimate the number of significant digits will be lost
            %Estimate that P-S=log(kappa) where P is input significant 
            %digits and S is number of reliable significant digits returned
            lostDigits=log10(kappa);
            fprintf('Likely to lose %.1f sig. digits\n\n', lostDigits)
        
            %Assign AFLAG based on how many significant digits will be lost
            %Matlab can store 16; 3 are desired in return; no more than
            %13 should be lost for success
            
            %AFLAG is stored as a property to access it later without 
            %running the entire analysis
            if lostDigits > 13
                self.AFLAG= 0;
            else
                self.AFLAG= 1;
            end
            AFLAG= self.AFLAG;
        end 
        
        
        %Method to compute the displacements and reactions of the structure
        function [DEFL, REACT]=ComputeDisplacementReactions(self)
            
            %Classify DOF
            [freeDOF, fixedDOF, knownDOF]= ClassifyDOF(self);
            
            %Obtaining the applied loads, displacements and fixed end
            %forces
            [Pf, Ps, Pn, FeFf, FeFs, FeFn, deltaN] = ...
                     CreateLoadVectors(self, freeDOF, fixedDOF, knownDOF);
            
            % Displacements at the Free DOFs
            deltaF= self.Kff\(Pf - FeFf - self.Kfn*deltaN);
            
            %Forces/Reactions at Fixed DOFs
            Rs= self.Ksf*deltaF + self.Ksn*deltaN + FeFs - Ps;
            
            %Forces/Reactions at Known DOFs
            Rn= self.Knf*deltaF + self.Knn*deltaN + FeFn - Pn;
            
            %Initializing DEFL to zeros of 6 rows (No. of DOFs in a node)
            %and total number of nodes as number of columns
            DEFL= zeros(6, self.nnodes);
            
            %Placing the values of DeltaF in the indeces corresponding to
            %the Free DOFs
            DEFL(freeDOF)= deltaF;
            
            %Placing the values of DeltaN in the indeces corresponding to
            %the Known DOFs
            DEFL(knownDOF)= deltaN;
            
            %Placing Zero as values in the indeces corresponding to
            %the Fixed DOFs as there will be zero displacements there
            DEFL(fixedDOF)= 0;
            
            %Transposing the DEFL matrix to have 6 columns (corresponding 
            %to 6 DOFs in a node) and number of nodes as the number of
            %rows- The format which MASTAN2 requires
            DEFL= DEFL';
            
            %Saving the Deflections as a property to access them later
            %without running the entire analysis
            self.DEFL= DEFL;
            
            %Initializing REACT to zeros of 6 rows (No. of DOFs in a node)
            %and total number of nodes as number of columns
            REACT= zeros(6, self.nnodes);
            
            %Placing Zero as values in the indeces corresponding to
            %the Free DOFs as there will be No reactions there
            REACT(freeDOF)= 0;
            
            %Placing the values of Pn in the indeces corresponding to
            %the Known DOFs
            REACT(knownDOF)= Rn;
            
            %Placing the values of Ps in the indeces corresponding to
            %the Fixed DOFs
            REACT(fixedDOF)= Rs;
            
            %Transposing the REACT matrix to have 6 columns (corresponding 
            %to 6 DOFs in a node) and number of nodes as the number of
            %rows- The format which MASTAN2 requires
            REACT= REACT';
            
            %Saving the Reactions as a property to access them later
            %without running the entire analysis
            self.REACT= REACT;

        end
        
        
        %Method to create the load vectors from the concentrated loads and
        %distributed loads and divide into pieces that correspond to 
        %degrees of freedom with free, fixed, and known displacements.
        function [Pf, Ps, Pn, FeFf, FeFs, FeFn, deltaN] = ...
                       CreateLoadVectors(self, freeDOF, fixedDOF, knownDOF)
            
            %Initializing the FeF matrix to all zeros
            FeF = zeros(self.nnodes*6,1);
            
            % Looping over all the element objects
            for i= 1:self.nele
                %Assembling the matrix of the fixed end forces of the 
                %structure by adding the values of the global FeF matrix of
                %the elements to the corresponding DOFs in the structure 
                %FeF matrix
                FeF(GetElementDOF(self.Elements(i)))=...
                    FeF(GetElementDOF(self.Elements(i)))...
                    + GetFixedEndForcesGlobal(self.Elements(i));
            end
            
            concen_t= self.concen';
            
            Pf= concen_t(freeDOF);
            Ps= concen_t(fixedDOF);
            Pn= concen_t(knownDOF);
            
            FeFf= FeF(freeDOF);
            FeFs= FeF(fixedDOF);
            FeFn= FeF(knownDOF);
            
            
            fixityTrans= self.fixity';
            deltaN= fixityTrans(knownDOF);

        end

        
        %Method to recover element forces from the element method and store
        %in the Mastan format
        function ELE_FOR=RecoverElementForces(self)

            %Initialize element forces to a matrix of zeroes
            ELE_FOR=zeros(self.nele, 12);

            %Loop over all elements, calling public method of element and
            %passing element object and deflections 
                for i=1:self.nele
                    ELE_FOR(i, 1:12)=ComputeForces(self.Elements(i), ...
                    [self.DEFL(self.ends(i,1), 1:6),self.DEFL(self.ends(i,2), 1:6)]);
                end
                
            %Saving the Element Forces as a property to access them later
            %without running the entire analysis
            self.ELE_FOR= ELE_FOR;
        end
    
    end
    
end 
    
