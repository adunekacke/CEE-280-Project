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
            Cw
            Zzz
            Zyy
            Ayy
            Azz
            E
            v
            Fy
            YldSurf
            Wt
            webdir
            DistribLoads
            Thermal
            Nodes
            Elements
            StructureStiffnessMatrix
            FeF
            P
    end
    
    % Public methods go here
    methods (Access = public)
        %% Constructor
        %  Replace XYZ by your initials before proceeding
        function self = ADSA_Analysis(nnodes, coord, concen, fixity, nele, ends,A,Izz,Iyy,J,Cw,Zzz,Zyy,...
		Ayy,Azz,E,v,Fy,YldSurf,Wt,webdir, w, thermal)
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
            self.Cw=Cw;
            self.Zzz=Zzz;
            self.Zyy=Zyy;
            self.Ayy=Ayy;
            self.Azz=Azz;
            self.E=E;
            self.v=v;
            self.Fy=Fy;
            self.YldSurf=YldSurf;
            self.Wt=Wt;
            self.webdir=webdir;
            self.DistribLoads=w;
            self.Thermal=thermal;
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
        function RunAnalysis(self)
        Kstruct = CreateStiffnessMatrix(self);
        [FeFstruct, Pstruct] = CreateLoadVectors(self);
        end
    end
    
    % Private methods go here
    methods (Access = private)
        
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
    end
end