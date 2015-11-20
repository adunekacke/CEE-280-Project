classdef ADSA_Element < handle

% Element class for a 3-dimensional framed structure

        %     element_nodes: A 2x1 vector of the Node Objects corresponding
        %                    to the start node and the end node
        %     A            : Element's cross sectional area
        %     Izz          : Element's moment of inertia about its local
        %                    z-z axis
        %     Iyy          : Element's moment of inertia about its local
        %                    y-y axis
        %     J            : Element's torsional constant
        %     Zzz          : Element's plastic section modulus about its
        %                    local z-z axis
        %     Zyy          : Element's plastic section modulus about its
        %                    local y-y axis
        %     Ayy          : Element's effective shear area along its local
        %                    y-y axis
        %     Azz          : Element's effective shear area along its local
        %                    z-z axis
        %     E            : Element's Young's Modulus
        %     v            : Element's material Poisson's ratio
        %     webDir       : A 3x1 vector of Element's unit web vector.  This is a unit vector
        %                          that defines the element's local y-y axis with respect
        %                          to the global coordinate system.  It is based on the
        %                          structure's undeformed geometry.
        %                              webDir(1) = x component of element's unit web vector
        %                              webDir(2) = y component of element's unit web vector
        %                              webDir(3) = z component of element's unit web vector
        %     w(i,1:3)         ==  element i's uniform load which references its
        %                            local coordinate system
        %                              w(i,1) = x component of uniform load
        %                              w(i,2) = y component of uniform load
        %                              w(i,3) = z component of uniform load
    
  %% Private properties go here
    properties (Access = private)
        elementNodes 
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
        distribLoad
        % The variables below are created within methods called by the constructor
        length
        gamma
        localStiffness
        globalStiffness
        elementDOF
        fixedEndForcesLocal
        fixedEndForcesGlobal
    end
    
    
  %% Public methods go here
    methods (Access = public)
        
        %% Constructor
        function self = ADSA_Element(elementNodes,A,Izz,Iyy,J,Zzz,Zyy,...
		                              Ayy,Azz,E,v,webDir,w)
            self.elementNodes=elementNodes;
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
            self.distribLoad=w;
            %Calling private method to compute and store the element length
            self.length=ComputeLength(self);  
            %Calling private method computes and store the element
            %transformation matrix
            self.gamma=ComputeTransformationMatrix(self); 
            %Calling the private method to compute and store the local and
            %global stiffness matrices of the element
            [self.localStiffness, self.globalStiffness]=ComputeElasticStiffnessMatrix(self); 
            %Storing the DOFs of the start and end nodes of the element in
            %a 12x1 column vector
            self.elementDOF=RetrieveDOF(self);
            %Computing and storing the local and global fixed end forces of
            %the element
            [self.fixedEndForcesLocal, self.fixedEndForcesGlobal]=ComputeFixedEndForces(self);
        end
        
        
        %% Getter functions
        
        %Get function for the global element stiffness matrix
        function globalStiffness = GetGlobalStiffness(self)
            globalStiffness = self.globalStiffness;
        end
        
        
        %Get function for the element DOFs
        function elementDOF = GetElementDOF(self)
            elementDOF = self.elementDOF;
        end
        
        
        %Get function for the global fixed end forces of the element
        function fixedEndForcesGlobal = GetFixedEndForcesGlobal(self)
            fixedEndForcesGlobal = self.fixedEndForcesGlobal;
        end
        
        
        %% Function to compute element forces in local coordinates
        function elementForces=ComputeForces(self, eleDelta)
            
            %Computing local forces using local stiffness, fixed end
            %forces, and gamma to transform global displacement at the
            %element's DOF's into element's local coordinates
            elementForces=(self.localStiffness*self.gamma*eleDelta')+...
                self.fixedEndForcesLocal;
        end
        
    end
    
    
  %% Private methods go here
    methods (Access = private)
        
        %% Compute the element's length
        function length=ComputeLength(self)
            
            %Call node method to get node coordinates at each end
            firstNode=GetNodeCoord(self.elementNodes(1)); 
            secondNode=GetNodeCoord(self.elementNodes(2)); 
            
            %Use pythagorean theorem to compute length of the element
            length=sqrt(((firstNode(1)-secondNode(1))^2)+((firstNode(2)-...
                      secondNode(2))^2)+((firstNode(3)-secondNode(3))^2));
                            
        end
        
        
        %% Compute the element's geometric transformation matrix
        function transformationMatrixGamma=ComputeTransformationMatrix(self)
            
            %Call node method to get node coordinates at each end
            firstNode=GetNodeCoord(self.elementNodes(1)); 
            secondNode=GetNodeCoord(self.elementNodes(2)); 
            
            %Create x' portion of small gamma matrix using vector projections
            xprime=(secondNode-firstNode)./self.length;
                
            %Create z' portion of small gamma using cross product of x' and
            %web direction unit vector
            zprime=cross(xprime,self.webDir);
            
            % Initializing the transformation matrix, Gamma to zeros
            transformationMatrixGamma=zeros(12);
            
            %Assembling the transformation matrix, Gamma
            transformationMatrixGamma(1:3,1:3)=[xprime;self.webDir;zprime];
            transformationMatrixGamma(4:6,4:6)=[xprime;self.webDir;zprime];
            transformationMatrixGamma(7:9,7:9)=[xprime;self.webDir;zprime];
            transformationMatrixGamma(10:12,10:12)=[xprime;self.webDir;zprime];
            
        end
        
        
        %% Compute the element's elastic stiffness matrix in local and global coordinates
        function [localStiffness, globalStiffness]=ComputeElasticStiffnessMatrix(self)
           
          %Store pieces of the local stiffness matrix (fewer computations
            %for efficiency)
            
            % Eta for Shear Deformations Stiffness
            etay=2*(1+self.v)*self.Izz/self.Ayy;  % Along Local Y-direction
            etaz=2*(1+self.v)*self.Iyy/self.Azz;  % Along Local Z-direction
            
            % Local Axial Deformation Stiffness
            oneone=self.A/self.length;
            
          % Stiffnesses due to Rotation and Bending Moments
            
            % Rotation about Y on one end due to Moment about Y on the
            % opposite end
            twoy=(self.Iyy/(self.length*...
                (((self.length^2/12))+etaz)))*(((self.length^2)/6)-etaz);
            
            % Rotation about Z on one end due to Moment about Z on the
            % opposite end
            twoz=(((self.Izz)/(self.length*...
                (((self.length^2/12))+etay)))*(((self.length^2)/6)-etay));
            
            % Rotation about Y on one end due to Moment about Y on the
            % same end
            foury=(self.Iyy/(self.length*...
                (((self.length^2/12))+etaz)))*(((self.length^2)/3)+etaz);
            
            % Rotation about Z on one end due to Moment about Z on the
            % same end
            fourz=(self.Izz/(self.length*...
                (((self.length^2/12))+etay))*(((self.length^2)/3)+etay));
            
            % Rotation about Y due to Force along Z OR Displacement along Z
            % due to Moment about Y
            sixy=((self.length/2)*...
                self.Iyy/(self.length*(((self.length^2/12))+etaz)));
            
            % Rotation about Z due to Force along Y OR Displacement along Y
            % due to Moment about Z
            sixz=((self.length/2)*...
                (self.Izz/(self.length*(((self.length^2/12))+etay))));
            
            % Displacement along Z due to Force along Z while causing
            % bending about Y
            twelvey=self.Iyy/(self.length*(((self.length^2/12))+etaz));
            
            % Displacement along Y due to Force along Y while causing
            % bending about Z
            twelvez=self.Izz/(self.length*(((self.length^2/12))+etay));
            
            % Torsional Stiffness about the Local X-axis
            torsion=self.J/(2*self.length*(1+self.v));

          %Create local member stiffness matrix
            
            localStiffness=(self.E)*[oneone     0      0     0     0       0    -oneone    0      0      0       0       0; ...
                                    0      twelvez   0     0     0     sixz      0    -twelvez  0      0       0     sixz; ...
                                    0         0   twelvey  0  -sixy      0       0       0  -twelvey   0     -sixy     0; ...
                                    0         0      0  torsion  0       0       0       0      0   -torsion   0       0; ...
                                    0         0    -sixy   0   foury     0       0       0     sixy    0      twoy     0; ...
                                    0        sixz    0     0     0     fourz     0    -sixz     0      0       0      twoz;...
                                    -oneone   0      0     0     0       0     oneone    0      0      0       0       0; ...
                                    0     -twelvez   0     0     0    -sixz      0    twelvez   0      0       0    -sixz; ...
                                    0         0  -twelvey  0   sixy      0       0       0   twelvey   0      sixy     0; ...
                                    0         0      0 -torsion  0       0       0       0      0    torsion   0       0; ...
                                    0         0    -sixy   0    twoy     0       0       0     sixy    0       foury   0; ...
                                    0        sixz    0     0     0      twoz     0     -sixz    0      0        0     fourz];
            
            %Tranform local stiffness matrix to global coordinates using
            %the transformation matrix.
            globalStiffness=transpose(self.gamma)*localStiffness*self.gamma;
            
        end
        
        
        %% Retrieving the element DOFs corresponding to the start and end
        %nodes of the element
        function elementDOF = RetrieveDOF(self)
            
            %The DOFs of the start and end nodes of the element are called
            %using the corresponding Get method of the Node class and then
            %the two sets of DOFs are appended in a single 12x1 column
            %vector
            elementDOF = [GetNodeDOF(self.elementNodes(1));...
                GetNodeDOF(self.elementNodes(2))];
            
        end
        
        
        %% Method to compute Fixed End Forces on the element due to
        %distributed loads
        function [fixedEndForcesLocal, fixedEndForcesGlobal]=ComputeFixedEndForces(self)
            
            %First, assemble a 12x1 matrix of local FEF's using general
            %equations.
            fixedEndForcesLocal=[-self.distribLoad(1)*self.length/2;
                                 -self.distribLoad(2)*self.length/2;
                                 -self.distribLoad(3)*self.length/2;
                                 0;
                                 (self.distribLoad(3)*self.length^2)/12;
                                 -(self.distribLoad(2)*self.length^2)/12;
                                 -self.distribLoad(1)*self.length/2;
                                 -self.distribLoad(2)*self.length/2;
                                 -self.distribLoad(3)*self.length/2;
                                 0;
                                 -(self.distribLoad(3)*self.length^2)/12;
                                 (self.distribLoad(2)*self.length^2)/12];
            
            %The FEF's in global coordinates are the transposed Gamma matrix times the local FEF's                 
            fixedEndForcesGlobal=self.gamma'*fixedEndForcesLocal;
            
        end 
        
    end
    
end
