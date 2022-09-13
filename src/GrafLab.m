%==========================================================================
%REFERENCES 
%==========================================================================
%
%Bucha, B., Janak, J., 2013. A MATLAB-based graphical user interface
%program for computing functionals of the geopotential up to ultra-high 
%degrees and orders. Computers and Geosciences 56, 186-196,
%http://dx.doi.org/10.1016/j.cageo.2013.03.012.
%
%Contact: blazej.bucha@stuba.sk
%
%
%
%==========================================================================
%USER MANUAL (see also the "docs" directory and the
%"https://github.com/blazej-bucha/graflab-cookbook" repository)
%==========================================================================
%
%This is a user manual for working with the GUI of the software. The manual 
%to run the software without the GUI is described later (look for the text 
%"function GrafLab").
%
%The GUI of the GrafLab is visually divided into three panels:
%
%--------------------------------------------------------------------------
%(i) GEOPOTENTIAL MODEL AND REFERENCE SYSTEM SELECTION: At
%first, the input GGM file or its error variance-covariance matrix
%must be imported using the "Browse. . ." button. The input
%GGM file must have one of the two standardized structures, see
%Table 1 and Table 2. In addition to the spherical harmonic
%coefficients, the input file may or may not contain the fifth and
%the sixth column with their standard deviations. GrafLab is also
%capable to read the input GGM file in a standard format defined by
%ICGEM (International Centre for Global Earth Models). In this case,
%the input file must have the suffix ".gfc" and the spherical harmonic
%coefficients sorted with respect to Table 1 or Table 2. If the input 
%GGM file has this particular structure, GrafLab uses the values of "GM" 
%and "R" from this file, and ignores the values of these two variables
%entered in the GUI. The input file may either be an ASCII file or a 
%binary MAT-file. In case of the GGM with a high maximum degree of SHE, 
%it is recommended to use the binary MAT-file, since it can be loaded 
%much faster.
%
%Table 1: Structure of the input GGM file - spherical harmonic 
%coefficients sorted primarily according to degrees.
%----------------------------------------
%  n   m       C_nm           S_nm
%----------------------------------------
%  2   0   -0.48417E-03    0.00000E+00
%  2   1   -0.20662E-09    0.13844E-08
%  2   2    0.24394E-05   -0.14003E-05
%  3   0    0.95716E-06    0.00000E+00
%----------------------------------------
%
%Table 2: Structure of the input GGM file - spherical harmonic 
%coefficients sorted primarily according to orders.
%----------------------------------------
%  n   m       C_nm           S_nm
%----------------------------------------
%  2   0   -0.48417E-03    0.00000E+00
%  3   0    0.95712E-06    0.00000E+00
%  4   0    0.53998E-06    0.00000E+00
%  5   0    0.68658E-07    0.00000E+00
%----------------------------------------
%
%The input ASCII file of error variance-covariance matrix
%must have the structure as shown in the Table 3. A binary
%MAT-file may be used to import error variance-covariance matrix
%as well. However, in this case, the empty arrays in Table 3 
%must be filled with zeroes or corresponding covariances.
%
%Table 3: Structure of the input file of error variance-covariance matrix - 
%spherical harmonic coefficients sorted primarily according to orders; 
%"n_min"=2; "n_max"=3; the column CS determines whether the variance and 
%covariances in the particular line are related to the coefficient "C_nm"
%(if "CS"=0) or to the coefficient "S_nm" (if "CS"=1).
%--------------------------------------------------------------------------------------------------------------------------
%  CS  n   m                variances and covariances of the spherical harmonic coefficients
%--------------------------------------------------------------------------------------------------------------------------
%  0   2   0   -4.31E-25
%  0   3   0   -2.11E-26-2.48E-25
%  0   2   1   -3.79E-28-1.15E-27-3.84E-25
%  1   2   1   -3.44E-28-4.67E-28-1.17E-27-4.16E-25
%  0   3   1   -1.99E-27-7.61E-29-2.98E-26-3.18E-28-2.48E-25
%  1   3   1   -1.44E-28-8.80E-29-3.42E-28-2.54E-26-3.16E-27-2.70E-25
%  0   2   2   -8.17E-27-1.72E-27-2.94E-28-3.67E-28-9.06E-29-1.08E-27 4.02E-25
%  1   2   2   -1.14E-27-2.94E-28-5.61E-29-3.86E-28-1.23E-27-1.50E-27-8.37E-28-4.25E-25
%  0   3   2   -9.38E-27-6.35E-27-1.08E-27-1.81E-27-7.12E-28-3.53E-28-3.30E-26-9.75E-29-3.07E-25
%  1   3   2   -1.27E-28-3.45E-27-1.59E-27-7.97E-28-1.75E-28-1.15E-28-5.51E-28-2.30E-26-2.78E-27-3.09E-25
%  0   3   3   -7.74E-28-1.36E-28-9.93E-27-5.50E-28-9.55E-28-3.25E-27-1.06E-27-8.60E-28-2.85E-29-1.58E-28-2.74E-25
%  1   3   3   -1.14E-27-2.19E-28-4.51E-28-1.26E-26-1.46E-28-4.90E-27-1.25E-28-1.76E-28-1.18E-28-5.22E-29-6.97E-29-2.74E-25
%--------------------------------------------------------------------------------------------------------------------------
%
%Most of GGMs have the same values of the geocentric gravitational
%constant and the radius of the reference sphere, therefore in this
%panel GrafLab automatically offers them for the computation.
%However, they may be simply replaced by the required values, 
%if necessary. Using the arrays "nmin" and "nmax", integer values
%in the intervals nmin 'in' <0,nmax> and nmax 'in' <2,M> may be
%entered (note that there are a few exceptions where the "nmin" value
%is fixed to 0 and cannot be changed, see Table 4 in (iii) "Calculated
%parameters and output selection" below). From the pop-up menu 
%"Ellipsoid", the normal gravity field generated by the 
%equipotential ellipsoid WGS84 (NIMA, 2000) or GRS80 (Moritz, 2000)
%can be selected.
%
%--------------------------------------------------------------------------
%(ii) POINT TYPE SELECTION: In the point type selection panel, it is
%possible to choose between the ellipsoidal or spherical type of the
%input coordinates. Next, one of three organizations of the evaluation 
%points must be specified by selecting the checkbox: "Grid", "Load data" 
%or "Point-wise". If the grid is selected, the minimum, maximum and 
%discretization step in the latitude (ellipsoidal or spherical) and 
%longitude directions must be entered. The array "Height above the 
%reference surface (m)" denotes the constant height of the grid above 
%the reference ellipsoid in the case of the ellipsoidal type of the 
%coordinates or above the reference sphere with the radius "R", defined
%by the GGM, in the case of the spherical coordinates. For the computation 
%on a regular grid, the lumped coeffcients approach is used. 
%
%To import the computational points from a data file, the "Load
%data" checkbox must be selected and subsequently, the data file
%must be imported using the Browse. . . button next to the checkbox.
%The data file may contain essentially an arbitrary number
%of lines and in every line of the file, the triplet of the 
%ellipsoidal/spherical coordinates (ellipsoidal/spherical latitude,
%longitude, both in degrees, and and ellipsoidal height/spherical radius
%in meters) must be given.
%
%After selecting the checkbox "Point-wise", an arbitrary point
%defined also by the triplet of the ellipsoidalspherical coordinates 
%can be entered manually using the arrays "Latitude (deg)", "Longitude (deg)", 
%"Ellipsoidal height/Spherical radius (m)". Type of the latitude in the
%array "Latitude (deg)" must correspond to the type of the input coordinates.
%In case of more points, the coordinates in each array must
%be separated by the comma or by the space. This point type
%selection is suitable if only a few points are to be determined,
%so there is no need to create a data file to import.
%
%In the last two mentioned cases of the point type selection,
%the lumped coeffcients approach cannot be applied due
%to irregular distribution of the points. Therefore we used two
%loops, one degree-depended and one order-depended.
%
%In each of the three above mentioned point type selections,
%the entries must be either in the form of floating point numbers
%with decimal dots or integer values, latitudes must be entered within the
%interval <-90 deg, 90 deg> and longitudes in the range <0 deg, 360 deg> or
%<-180 deg, 180 deg>.
%
%--------------------------------------------------------------------------
%(iii) CALCULATED PARAMETERS AND OUTPUT SELECTION: Using the
%four pop-up menus on the left side of this panel, user can simply
%choose, which functionals of the geopotential are to be computed.
%Note that at least one and maximum four functionals
%may be computed simultaneously. The summary of the functionals that can
%be computed in GrafLab is shown in Table 4. In order to stay brief, we do
%not introduce here the mathematical formulae for computing each
%functional, but these can be found in the pdf file
%"Definition_of_functionals_of_the_geopotential_used_in_GrafLab_software.pdf". 
%For evaluating disturbing and gravitational tensor in the LNOF, we
%used the non-singular expressions, which can be found e.g. in Petrovskaya
%and Vershkov (2006). For practical reasons, we slightly modified these
%formulae. The modified formulae can be found in the same pdf file.
%
%Table 4: Functionals of the geopotential available in GrafLab.
%Explanation of the symbols in the table: "V" - gravitational potential,
%"W" - gravity potential, "gravity_vector_gX_gY_gZ" - gravity vector, 
%"T" - disturbing potential, "delta g" - gravity disturbance, 
%"DELTA g" - gravity anomaly, "xi" - north-south component of deflection of the vertical, 
%"eta" - east-west component of deflection of the vertical, 
%"THETA" - total deflection of the vertical, "N" - geoid undulation, 
%"zeta_Ell" - generalized height anomaly, "zeta" - height anomaly; the subscript "sa" denotes the spherical
%approximation of the functional; ("r", "theta", "lambda") stands for the
%spherical coordinates; ("x","y","z") denotes the coordinates in the local
%north-oriented reference frame; the subscripts "r", "theta", "lambda", 
%"x", "y", "z" and their combinations stand for the derivatives of the 
%functionals with respect to the particular coordinate; the number in the 
%superscript denotes computational demand (computation time of the 
%functional and memory usage during the computation) - (1) small, 
%(2) medium, (3) high; (*) denote the functionals for which the 
%value of nmin cannot be larger than 0.
%-------------------------------------------------------------------------------------------------------------------
%          Actual field                    Disturbing field          Geometrical characteristics of the actual field
%-------------------------------------------------------------------------------------------------------------------
%              V(1)                             T(1)                           xi(2)
%   V_rr(3),V_phiphi(3),V_ll(3)     T_rr(3),T_phiphi(3),T_ll(3)                eta(1)
%   V_rphi(3),V_rl(3),V_phil(3)     T_rphi(3),T_rl(3),T_phil(3)                THETA(2)
%   V_xx(2),V_yy(2),V_zz(2)         T_xx(2),T_yy(2),T_zz(2)                    (*)N(2)
%(*)V_xy(3),(*)V_xz(3),(*)V_yz(3)   (*)T_xy(3),(*)T_xz(3),(*)T_yz(3)           zeta_Ell(1)
%              W(1)                          (*)delta g(2)                     (*)zeta(2)
%   gravity_vector_gX_gY_gZ(2)               delta g_sa(1)
%              g_sa(1)                       DELTA g_sa(1)
%              W_rr(1)                          T_rr(1)
%-------------------------------------------------------------------------------------------------------------------
%
%To compute geoid undulation "N" and height anomaly "zeta", the
%digital terrain model, e.g. DTM2006.0 (Pavlis et al., 2007),
%must be imported. Only one particular structure of the DTM
%file, shown in Table 1, can be recognized by the GrafLab.
%If these two functionals are to be computed, immediately after
%clicking the "OK" button, the dialog window from which the input
%DTM file must be imported will appear.
%
%Each functional of the geopotential may be evaluated using
%any of the three approaches for computing fnALFs except for the 
%gravitational and disturbing tensors in the LNOF. Since these 
%non-singular expressions have been slightly modified, the modified 
%forward column method combined with Horner's scheme is not effcient 
%for the new formulae and therefore it was not used in this case.
%
%In order to evaluate the commission errors of the functionals,
%the "Commission error" check box must be selected. GrafLab allows
%to compute the commission errors of the each above mentioned functional 
%except for the gravitational and disturbing tensors in the LNOF: 
%"T_xx"; "T_yy"; "T_zz"; "T_xy"; "T_xz"; "T_yz"; "V_xx"; "V_yy"; "V_zz"; 
%"V_xy"; "V_xz"; "V_yz". One should keep in mind that the evaluation of 
%commission error has much higher requirements on the PC, because of the 
%large size of the error variance-covariance matrix. In general it means 
%that the maximum degree "M" is reduced from thousands and hundreds to tens, 
%and number of computing points have to be decreased as well.
%
%By clicking the button "Computation of fnALFs", a new dialog
%window will appear, in which user may choose one of the three 
%approaches for evaluating values of fnALFs.
%
%If the Mapping toolbox of MATLAB is installed, computed data may be 
%depicted on a map using automatically selected cartographic projection 
%(e.g. pseudocylindical Robinson projection, equidistant conic projection, 
%equidistant azimuthal projection, ...). By clicking the button "Display 
%data settings", another dialog window will appear. Here, user can set 
%up the required output parameters of the exported map. This option is 
%available only if the computation on a regular grid has been chosen.
%
%The button "Output folder and file" permits to specify the output
%folder and prefix of the all exported files, i.e. without any
%suffix (e.g. "Prefix"). The data file (e.g. "Prefix.txt") with the computed
%data may be created by selecting the checkbox "Export data". 
%The report file, which contains the informations about the
%computation, may be created by selecting the
%"Export report" checkbox. This file automatically obtains name
%with the suffix "_Report.txt", e.g. "Prefix_Report.txt". If the "Display
%data" checkbox has been selected, GrafLab creates also a
%graphical file (or files, depending on the number of computing
%functionals) according to chosen graphic file format (bmp, emf,
%eps, jpeg, pdf, png or tiff).
%
%When all the required input parameters and input files have
%been entered, after clicking the "OK" button, the computation
%will start. On the left from this button, there is a status line,
%which provides short explanations during the whole computational
%process ("Loading GGM file..." , current value of the variable
%m in the order-dependent loop, "Displaying data..." , ... ), so
%one can clearly see in which part of the computation is GrafLab.
%After successful computation, the status "Computation has been
%finished" will appear. If any of the input parameters or input
%files have been entered in a wrong format, GrafLab will open a
%message dialog or error dialog with description of the error.

function output=GrafLab(vstpar,... %To work without the GUI, set this variable to the string 'OK'
    ...
    ...                  %GLOBAL GEOPOTENTIAL MODEL AND REFERENCE SYSTEM SELECTION
    ...                  %==========================================================================
    GM,...               %Geocentric gravitational constant of the GGM (m^3*s^-2)
    R,...                %Radius of the reference sphere of the GGM (m)
    nmin,...             %Minimum degree of the spherical harmonic expansion
    nmax,...             %Maximum degree of the spherical harmonic expansion
    ...                  %(to automatically use the maximum degree of the imported
    ...                  %global geopotential model, set this variable to the string 'nmaxGGM')
    ellipsoid,...        %Reference ellipsoid: this variable can either be
    ...                  %a scalar (1 - GRS80, 2 - WGS84) or a vector with
    ...                  %5 elements [GM a e C20 omega], where "GM" is the 
    ...                  %geocentric gravitational constant (m^3*s^-2), "a"
    ...                  %is the semimajor axis (m), "e" is the first eccentricity,
    ...                  %"C20" is the fully normalized C20 and "omega"
    ...                  %is the angular velocity (rad*sec^-1)
    GGM_path,...         %Path (absolute or relative) to the GGM file (or to the covariance matrix file if "Functional_or_Commission==1", see below)
    ...                  %==========================================================================
    ...
    ...                  %POINT TYPE SELECTION
    ...                  %==========================================================================
    coord,...            %Type of input coordinates: 0 - ellipsoidal, 1 - spherical
    point_type,...       %Point type: 0 - Grid, 1 - Load data, 2 - Point-wise
    ...
    ...                  %Specifications for the GRID MODE (i.e., if "point_type==0"; otherwise, the variables
    ...                  %"lat_min", ..., "h" are ignored and can be set to, e.g., "[]").
    ...                  %--------------------------------------------------------------------------
    ...                  %The grid can be defined in two ways, here denoted as
    ...                  %i) and ii). A mix of i) and ii) is not allowed.
    lat_min,...          %i)  Latitude minimum in degrees (scalar)
    ...                  %ii) Vector defining parallels of the grid in degrees
    lat_step,...         %i)  Latitude step in degrees (scalar)
    ...                  %ii) String 'empty'
    lat_max,...          %i)  Latitude maximum in degrees (scalar)
    ...                  %ii) String 'empty'
    lon_min,...          %i)  Longitude minimum in degrees (scalar)
    ...                  %ii) Vector defining meridians of the grid in degrees
    lon_step,...         %i)  Longitude step in degrees (scalar)
    ...                  %ii) String 'empty'
    lon_max,...          %i)  Longitude maximum in degrees (scalar)
    ...                  %ii) String 'empty'
    h,...                %Height above the reference surface in metres (scalar).
    ...                  %In case of ellipsoidal input coordinates ("coord==0"), the reference ellipsoid (e.g., GRS80
    ...                  %or WGS84) plays the role of the reference surface, and therefore "h" is the ellipsoidal
    ...                  %height.
    ...                  %In case of spherical input coordinates ("coord==1"), the reference surface is a sphere
    ...                  %with the radius "R" (see above), and therefore "h" represents the height above this sphere
    ...                  %measured in the radial direction.
    ...                  %--------------------------------------------------------------------------
    ...
    ...                  %Specifications for the LOAD DATA MODE (i.e., if "point_type==1"; otherwise, the variable
    ...                  %"Input_data_path" is ignored and can be set to, e.g., "[]").
    ...                  %--------------------------------------------------------------------------
    Input_data_path,...  %Path (absolute or relative) to the file specifying computation points.
    ...                  %--------------------------------------------------------------------------
    ...
    ...                  %Specifications for the POINT-WISE MODE (i.e., if "point_type==2"; otherwise, the variables
    ...                  %"lat", "lon" and "h2" are ignored and can be set to, e.g., "[]").
    ...                  %--------------------------------------------------------------------------
    lat,...              %Latitudes in degrees (column vector or scalar)
    lon,...              %Longitues in degrees (column vector or scalar)
    h2,...               %Ellipsoidal heights (if "coord==0") or spherical radii (if "coord==1"). Both in metres
    ...                  %--------------------------------------------------------------------------
    ...                  %==========================================================================
    ...
    ...                  %CALCULATED PARAMETERS AND OUTPUT SELECTION
    ...                  %==========================================================================
    Output_path,...      %Path (absolute or relative) to the output data file (without any suffix of the output file)
    Functional_or_Commission,... %Computation of: 0 - functional, 1 - commission error 
    Functional,...       %Column vector defining functionals to be computed (see the list below). The total number of elements of this vector cannot exceed 4.
    fnALFs,...           %Computation of fnALFs using: 1 - standard forward column method, 2 - modified forward column method, 3 - extended-range arithmetic
    DTM_path,...         %Path (absolute or relative) to the spherical harmonic coefficients of a DTM (if "Geoid undulation" or "Height anomaly" is to be computed, a path to the file with spherical harmonic coefficients of a digital terrain model has to be specified)
    ...
    Export_data_txt,...  %Export data into a txt file: 0 - No, 1 - Yes
    Export_report,...    %Export report:               0 - No, 1 - Yes
    Export_data_mat,...  %Export data into a mat file: 0 - No, 1 - Yes
    ...
    ...                  %Display data settings
    ...                  %--------------------------------------------------------------------------
    Display_data,...     %Display data (possible only if "point_type==0"): 
    ...                  %0 - No
    ...                  %1 - Yes, uses Mapping Toolbox (slow, but the result is a nice countour map
    ...                  %    in a cartographic projection)
    ...                  %2 - Yes, uses the "imagesc" function (fast, 
    ...                  %    but the result is plotted only as a matrix)
    ...                  %If "Display_data==0" or "isempty(Display_data)",
    ...                  %the variables "Graphic_format", ..., "DPI" are ignored
    ...                  %and can be set to, e.g., "[]"
    Graphic_format,...   %Graphic file format: 1 - bmp, 2 - emf, 3 - eps, 4 - jpeg, 5 - pdf, 6 - png, 7 - tiff
    Colormap,...         %Colormap: 1 - jet, 2 - hsv, 3 - hot, 4 - cool, 5 - spring, 6 - summer, 7 - autumn, 8 - winter, 9 - gray, 10 - bone, 11 - copper, 12 - pink, 13 - lines
    Number_of_colors,... %Number of colors of the selected colormap
    DPI,...              %Hopefully, this variable is self-explanatory :)
    ...                  %--------------------------------------------------------------------------
    ...                  %==========================================================================
    Status_bar)          %Display the computational progress in the command window: 0 - No, 1 - Yes
%
%The list of the functionals that can be computed (the numbers assigned to
%the functionals are the input to the variable "Functional"; see above):
%2  - Deflection of the vertical eta (east-west component)
%3  - Deflection of the vertical xi (north-south component)
%4  - Deflection of the vertical Theta (total)
%5  - Disturbing potential
%6  - Disturbing tensor (T_rr,   T_phiphi,  T_lambdalambda)
%7  - Disturbing tensor (T_rphi, T_rlambda, T_philambda)
%8  - Disturbing tensor in the LNOF (Txx, Tyy, Tzz)
%9  - Disturbing tensor in the LNOF (Txy, Txz, Tyz)
%10 - Geoid undulation
%11 - Gravitational potential
%12 - Gravitational tensor (V_rr,   V_phiphi,  V_lambdalambda)
%13 - Gravitational tensor (V_rphi, V_rlambda, V_philambda)
%14 - Gravitational tensor in the LNOF (V_xx, V_yy, V_zz)
%15 - Gravitational tensor in the LNOF (V_xy, V_xz, V_yz)
%16 - Gravity vector in the LNOF (g_X, g_Y, g_Z)
%17 - Gravity in spherical approximation (magnitude of 
%     the gravity vector in spherical approximation)
%18 - Gravity potential
%19 - Gravity anomaly in spherical approximation
%20 - Gravity disturbance (difference between magnitudes
%     of the gravity vector and of the normal gravity vector)
%21 - Gravity disturbance in spherical approximation
%22 - Height anomaly ell (approximation of height anomaly or geoid)
%23 - Height anomaly
%24 - Second radial derivative of disturbing potential
%25 - Second radial derivative of gravity potential

R1=0.8; G1=0.8; B1=0.8;
R2=0.95; G2=0.95; B2=0.95;

font='Helvetica';
screenres=get(0,'screensize');

if nargin==0  
     
    %Main window 
    M=figure('units','pixels','numbertitle','off','name','GrafLab 2.1.4',...
        'color',[R1 G1 B1],'position',[screenres(3)/2-700/2 screenres(4)/2-640/2 700 640],...
        'tag','okno','menubar','none');      
    a=0.01; b=0.04; c=0.045; d=-0.02;
    
    %Panels
    %======================================================================
    %Geopotential model and reference system selection panel
	uipanel('Units','normalized','position',[0.03 0.77 0.935 0.21],...
        'backgroundcolor',[R1 G1 B1],'HighlightColor',[R2 G2 B2],'tag',...
        'GMaRSSpanel');
    
    %Point type selection panel
	uipanel('Units','normalized','position',[0.03 0.355 0.935 0.4],...
        'backgroundcolor',[R1 G1 B1],'HighlightColor',[R2 G2 B2],'tag',...
        'PTSpanel');
    
    %Calculated parameters and output selection panel
	uipanel('Units','normalized','position',[0.03 0.08 0.935 0.26],...
        'backgroundcolor',[R1 G1 B1],'HighlightColor',[R2 G2 B2],'tag',...
        'CPaOSpanel');
    
    %Geopotential model and reference system selection
    %======================================================================
    uicontrol('Units','normalized','position',[0.05 0.865+a 0.14 0.035],...
        'style','pushbutton','string','Browse...','tag','GGM',...
        'callback','GrafLab import_GGM'); %Browse... button
    uicontrol('Units','normalized','position',[0.605 0.865+a 0.32 0.035],...
        'style','checkbox','string','Use maximum degree of GGM',...
        'value',1,'backgroundcolor',[R1 G1 B1],'tag','use'); %Use maximum degree of GGM
    uicontrol('Units','normalized','position',[0.05 0.91+a 0.375 0.025],...
        'style','text','string','Global geopotential model of the Earth',...
        'backgroundcolor',[R1 G1 B1]); %Text Global geopotential model of the Earth
    uicontrol('Units','normalized','position',[0.05 0.82+a 0.25 0.025],...
        'style','text','string','GM of GGM (m3*s-2)',...
        'backgroundcolor',[R1 G1 B1]); %Text GM of GGM (m3*s-2)
    uicontrol('Units','normalized','position',[0.05 0.78+a 0.25 0.035],...
        'style','edit','string','3986004.415E+8','backgroundcolor',...
        [R2 G2 B2],'tag','GM'); %Value of GM
    uicontrol('Units','normalized','position',[0.33 0.82+a 0.25 0.025],...
        'style','text','string','R of GGM (m)',...
        'backgroundcolor',[R1 G1 B1],'tag','R_text'); %Text R of GGM (m)
    uicontrol('Units','normalized','position',[0.33 0.78+a 0.25 0.035],...
        'style','edit','string','6378136.3','backgroundcolor',[R2 G2 B2],...
        'tag','R'); %Value of R
    uicontrol('Units','normalized','position',[0.605 0.82+a 0.08 0.025],...
        'style','text','string','nmin','backgroundcolor',[R1 G1 B1],...
        'tag','text_nmin'); %Text nmin
    uicontrol('Units','normalized','position',[0.605 0.78+a 0.08 0.035],...
        'style','edit','string','0','backgroundcolor',[R2 G2 B2],...
        'tag','nmin');  %Value of nmin
    uicontrol('Units','normalized','position',[0.71 0.82+a 0.08 0.025],...
        'style','text','string','nmax','backgroundcolor',[R1 G1 B1],...
        'tag','text_nmax'); %Text nmax
    uicontrol('Units','normalized','position',[0.71 0.78+a 0.08 0.035],...
        'style','edit','backgroundcolor',[R2 G2 B2],'tag','nmax'); %Value of nmax
    uicontrol('Units','normalized','position',[0.81 0.82+a 0.13 0.025],...
        'style','text','string','Ellipsoid','backgroundcolor',[R1 G1 B1],...
        'tag','ell_text'); %Text Ellipsoid
    uicontrol('Units','normalized','position',[0.815 0.717+a 0.13 0.1],...
        'style','popup','string','GRS80|WGS84','backgroundcolor',...
        [R2 G2 B2],'tag','ell'); %Ellipsoid - pop-up menu
    uicontrol('units','normalized','position',[0.22 0.865+a 0.36 0.035],...
        'style','edit','tag','nameGGM','enable','off'); %Name of the imported GGM file
    uicontrol('Units','normalized','position',[0.04 0.03 0.3 0.025],...
        'style','text','backgroundcolor',[R1 G1 B1],'tag','hlasky'); %Status line

    %Point type seleection
    %======================================================================
    
    %Text Input coordinates
    uicontrol('Units','normalized','position',[0.05 0.66+b 0.28 0.025],...
        'style','text','string','Type of the input coordinates:','backgroundcolor',...
        [R1 G1 B1]);
    
    %Radio button group (Ellipsoida, Spherical)
    c0 = uibuttongroup('visible','on','units','normalized',...
        'Position',[0.37 0.628+b 0.4 0.06],'bordertype','none',...
        'backgroundcolor',[R1 G1 B1],'tag','coordinates'); 
    c1 = uicontrol('units','normalized','Style','Radio','pos',...
        [0.0 0.5 0.4 0.5],'parent',c0,'HandleVisibility','on',...
        'backgroundcolor',[R1 G1 B1],'tag','rbutton1coord'); %Radio button: Ellipsoidal
    set(c1,'String','Ellipsoidal','fontname',font,'fontsize',10);
    c2 = uicontrol('units','normalized','Style','Radio','pos',...
        [0.388 0.5 0.4 0.5],'parent',c0,'HandleVisibility','on',...
        'backgroundcolor',[R1 G1 B1],'tag','rbutton2coord'); %Radio button: Spherical
    set(c2,'String','Spherical','fontname',font,'fontsize',10); 
    
    %Radio button group (Grid, Load data, Point-wise)
    d0 = uibuttongroup('visible','on','units','normalized',...
        'Position',[0.05 0.534+b 0.7 0.06],'bordertype','none',...
        'backgroundcolor',[R1 G1 B1],'tag','pointtype');
    d1 = uicontrol('units','normalized','Style','Radio','pos',...
        [0.0 0.5 0.4 0.5],'parent',d0,'HandleVisibility','on',...
        'backgroundcolor',[R1 G1 B1],'tag','rbutton1type'); %Radio button: Grid
    set(d1,'String','Grid','fontname',font,'fontsize',10);
    d2 = uicontrol('units','normalized','Style','Radio','pos',...
        [0.23 0.5 0.4 0.5],'parent',d0,'HandleVisibility','on',...
        'backgroundcolor',[R1 G1 B1],'tag','rbutton2type'); %Radio button: Load data
    set(d2,'String','Load data','fontname',font,'fontsize',10);
    d3 = uicontrol('units','normalized','Style','Radio','pos',...
        [0.735 0.5 0.4 0.5],'parent',d0,'HandleVisibility','on',...
        'backgroundcolor',[R1 G1 B1],'tag','rbutton3type'); %Radio button: Point-wise
    set(d3,'String','Point-wise','fontname',font,'fontsize',10);
    
    uicontrol('Units','normalized','position',[0.37 0.561+b 0.145 0.035],...
        'style','pushbutton','string','Browse...','tag','import',...
        'callback','GrafLab input'); %Browse... button
    
    uicontrol('units','normalized','position',[0.21 0.64+a 0.305 0.035],...
        'style','edit','tag','LoadFile','enable','off'); %Name of the input file with computation points

    %Grid
    %----------------------------------------------------------------------
    
    %phi
    uicontrol('Units','normalized','position',[0.05 0.48+b 0.145 0.035],...
        'style','edit','backgroundcolor',[R2 G2 B2],'tag','fimin'); %phi min
    uicontrol('Units','normalized','position',[0.21 0.48+b 0.145 0.035],...
        'style','edit','backgroundcolor',[R2 G2 B2],'tag','fistep'); %phi step
    uicontrol('Units','normalized','position',[0.37 0.48+b 0.145 0.035],...
        'style','edit','backgroundcolor',[R2 G2 B2],'tag','fimax'); %phi max

    %lambda
    uicontrol('Units','normalized','position',[0.05 0.407+b 0.145 0.035],...
        'style','edit','backgroundcolor',[R2 G2 B2],'tag','lambdamin'); %lambda min
    uicontrol('Units','normalized','position',[0.21 0.407+b 0.145 0.035],...
        'style','edit','backgroundcolor',[R2 G2 B2],'tag','lambdastep'); %lambda step
    uicontrol('Units','normalized','position',[0.37 0.407+b 0.145 0.035],...
        'style','edit','backgroundcolor',[R2 G2 B2],'tag','lambdamax'); %lambda max

    %h
    uicontrol('Units','normalized','position',[0.05 0.335+b 0.465 0.035],...
        'style','edit','backgroundcolor',[R2 G2 B2],'tag','hgrid'); %h
    
    %Text for latitudes
    uicontrol('Units','normalized','position',[0.05 0.52+b 0.145 0.025],...
        'style','text','string','Lat. min (deg)','backgroundcolor',...
        [R1 G1 B1],'tag','fimin_string'); %Text Lat. min (deg)
    uicontrol('Units','normalized','position',[0.21 0.52+b 0.145 0.025],...
        'style','text','string','Lat. step (deg)','backgroundcolor',...
        [R1 G1 B1],'tag','fistep_string'); %Text Lat. step (deg)
    uicontrol('Units','normalized','position',[0.37 0.52+b 0.145 0.025],...
        'style','text','string','Lat. max (deg)','backgroundcolor',...
        [R1 G1 B1],'tag','fimax_string'); %Text Lat. max (deg)
    
    %Text for longitudes
    uicontrol('Units','normalized','position',[0.05 0.445+b 0.145 0.025],...
        'style','text','string','Lon. min (deg)','backgroundcolor',...
        [R1 G1 B1]); %Text Lon. min (deg)
    uicontrol('Units','normalized','position',[0.21 0.445+b 0.145 0.025],...
        'style','text','string','Lon. step (deg)','backgroundcolor',...
        [R1 G1 B1]); %Text Lon. step (deg)
    uicontrol('Units','normalized','position',[0.37 0.445+b 0.145 0.025],...
        'style','text','string','Lon. max (deg)','backgroundcolor',...
        [R1 G1 B1],'tag','diskcheck'); %Text Lon. max (deg)
        
    %Text Height above the reference surface (m)
    uicontrol('Units','normalized','position',[0.05 0.37+b 0.465 0.025],...
        'style','text','string','Height above the reference surface (m)',...
        'backgroundcolor',[R1 G1 B1],'tag','h_string'); 

    %Point-wise
    %----------------------------------------------------------------------
    
    uicontrol('Units','normalized','position',[0.585+d 0.49+b 0.38 0.035],...
        'style','edit','backgroundcolor',[R2 G2 B2],'tag','fi'); %phi
    uicontrol('Units','normalized','position',[0.585+d 0.413+b 0.38 0.035],...
        'style','edit','backgroundcolor',[R2 G2 B2],'tag','lambda'); %lambda
    uicontrol('Units','normalized','position',[0.585+d 0.335+b 0.38 0.035],...
        'style','edit','backgroundcolor',[R2 G2 B2],'tag','hdisk'); %h

    uicontrol('Units','normalized','position',[0.585+d 0.53+b 0.38 0.025],...
        'style','text','string','Latitude (deg)','backgroundcolor',[R1 G1 B1]); %Text Latitude (deg)
    uicontrol('Units','normalized','position',[0.585+d 0.45+b 0.38 0.025],...
        'style','text','string','Longitude (deg)','backgroundcolor',[R1 G1 B1]); %Text Longitude deg)
    uicontrol('Units','normalized','position',[0.585+d 0.375+b 0.38 0.025],...
        'style','text','string','Ellipsoidal height/Spherical radius (m)',...
        'backgroundcolor',[R1 G1 B1]); %Ellipsoidal height/Spherical radius (m)

    %Calculated parameters and output selection
    %======================================================================
    
    %The first functional
    uicontrol('Units','normalized','position',[0.05 0.16+c 0.385 0.1],...
        'style','popup','string',...
        '|Deflection_of_the_vertical_eta|Deflection_of_the_vertical_xi|Deflection_of_the_vertical_Theta|Disturbing_potential|Disturbing_tensor_Trr_Tpp_Tll|Disturbing_tensor_Trp_Trl_Tpl|Disturbing_tensor_Txx_Tyy_Tzz|Disturbing_tensor_Txy_Txz_Tyz|Geoid_undulation|Gravitational_potential|Gravitational_tensor_Vrr_Vpp_Vll|Gravitational_tensor_Vrp_Vrl_Vpl|Gravitational_tensor_Vxx_Vyy_Vzz|Gravitational_tensor_Vxy_Vxz_Vyz|Gravity_vector_gX_gY_gZ|Gravity_sa|Gravity_potential|Gravity_anomaly_sa|Gravity_disturbance|Gravity_disturbance_sa|Height_anomaly_ell|Height_anomaly|2nd_radial_der_of_disturbing_potential|2nd_radial_der_of_gravity_potential',...
        'backgroundcolor',[R2 G2 B2],'tag','P1');
    %The second functional
	uicontrol('Units','normalized','position',[0.05 0.105+c 0.385 0.1],...
        'style','popup','string',...
        '|Deflection_of_the_vertical_eta|Deflection_of_the_vertical_xi|Deflection_of_the_vertical_Theta|Disturbing_potential|Disturbing_tensor_Trr_Tpp_Tll|Disturbing_tensor_Trp_Trl_Tpl|Disturbing_tensor_Txx_Tyy_Tzz|Disturbing_tensor_Txy_Txz_Tyz|Geoid_undulation|Gravitational_potential|Gravitational_tensor_Vrr_Vpp_Vll|Gravitational_tensor_Vrp_Vrl_Vpl|Gravitational_tensor_Vxx_Vyy_Vzz|Gravitational_tensor_Vxy_Vxz_Vyz|Gravity_vector_gX_gY_gZ|Gravity_sa|Gravity_potential|Gravity_anomaly_sa|Gravity_disturbance|Gravity_disturbance_sa|Height_anomaly_ell|Height_anomaly|2nd_radial_der_of_disturbing_potential|2nd_radial_der_of_gravity_potential',...
        'backgroundcolor',[R2 G2 B2],'tag','P2');    
    %The third functional
	uicontrol('Units','normalized','position',[0.05 0.05+c 0.385 0.1],...
        'style','popup','string',...
        '|Deflection_of_the_vertical_eta|Deflection_of_the_vertical_xi|Deflection_of_the_vertical_Theta|Disturbing_potential|Disturbing_tensor_Trr_Tpp_Tll|Disturbing_tensor_Trp_Trl_Tpl|Disturbing_tensor_Txx_Tyy_Tzz|Disturbing_tensor_Txy_Txz_Tyz|Geoid_undulation|Gravitational_potential|Gravitational_tensor_Vrr_Vpp_Vll|Gravitational_tensor_Vrp_Vrl_Vpl|Gravitational_tensor_Vxx_Vyy_Vzz|Gravitational_tensor_Vxy_Vxz_Vyz|Gravity_vector_gX_gY_gZ|Gravity_sa|Gravity_potential|Gravity_anomaly_sa|Gravity_disturbance|Gravity_disturbance_sa|Height_anomaly_ell|Height_anomaly|2nd_radial_der_of_disturbing_potential|2nd_radial_der_of_gravity_potential',...
        'backgroundcolor',[R2 G2 B2],'tag','P3');   
    %The fourth functional
	uicontrol('Units','normalized','position',[0.05 -0.005+c 0.385 0.1],...
        'style','popup','string',...
        '|Deflection_of_the_vertical_eta|Deflection_of_the_vertical_xi|Deflection_of_the_vertical_Theta|Disturbing_potential|Disturbing_tensor_Trr_Tpp_Tll|Disturbing_tensor_Trp_Trl_Tpl|Disturbing_tensor_Txx_Tyy_Tzz|Disturbing_tensor_Txy_Txz_Tyz|Geoid_undulation|Gravitational_potential|Gravitational_tensor_Vrr_Vpp_Vll|Gravitational_tensor_Vrp_Vrl_Vpl|Gravitational_tensor_Vxx_Vyy_Vzz|Gravitational_tensor_Vxy_Vxz_Vyz|Gravity_vector_gX_gY_gZ|Gravity_sa|Gravity_potential|Gravity_anomaly_sa|Gravity_disturbance|Gravity_disturbance_sa|Height_anomaly_ell|Height_anomaly|2nd_radial_der_of_disturbing_potential|2nd_radial_der_of_gravity_potential',...
        'backgroundcolor',[R2 G2 B2],'tag','P4');

    %----------------------------------------------------------------------
     
    %Commission error checkbox
    uicontrol('Units','normalized','position',[0.46 0.228+c 0.2 0.025],...
        'style','checkbox','string','Commission error','backgroundcolor',[R1 G1 B1],...
        'tag','STD'); %STD  
    
    %Computation of fnALFs
    uicontrol('Units','normalized','position',[0.46 0.167+c 0.22 0.038],...
        'style','pushbutton','string','Computation of fnALFs',...
        'callback','GrafLab fnALFs','tag','fnALFs'); 
    
    %Display data settings
    uicontrol('Units','normalized','position',[0.46 0.112+c 0.22 0.038],...
        'style','pushbutton','string','Display data settings',...
        'callback','GrafLab Display_data_settings','tag','DDS');
    
    %Output folder and file
    uicontrol('Units','normalized','position',[0.46 0.057+c 0.22 0.038],...
        'style','pushbutton','string','Output folder and file',...
        'callback','GrafLab Output_folder','tag','outfolder');
        
    %Export data checkbox
    uicontrol('Units','normalized','position',[0.725+d 0.228+c 0.24 0.025],...
        'style','checkbox','string','Export data','tag','export',...
        'callback','GrafLab output','backgroundcolor',[R1 G1 B1],'value',1);
    
    %Export report checkbox
    uicontrol('Units','normalized','position',[0.725+d 0.175+c 0.24 0.025],...
        'style','checkbox','string','Export report','backgroundcolor',...
        [R1 G1 B1],'tag','report','value',1);
    
    %Export data in *.mat
    uicontrol('Units','normalized','position',[0.725+d 0.12+c 0.24 0.025],...
        'style','checkbox','string','Export data in *.mat','backgroundcolor',...
        [R1 G1 B1],'tag','datamat','value',0);
    
    %Name of the output file
    uicontrol('units','normalized','position',[0.725+d 0.057+c 0.24 0.035],...
        'style','edit','tag','OutFile','enable','off');
    
    %OK button
    uicontrol('Units','normalized','position',[0.35 0.015 0.13 0.05],...
        'style','pushbutton','string','OK','Callback','GrafLab OK');
    
    %Close button
	uicontrol('Units','normalized','position',[0.55 0.015 0.13 0.05],...
        'style','pushbutton','string','Close','Callback','GrafLab Close');

    %Set font to Cambria
    set(get(M,'children'),'fontname',font,'fontsize',10);  
    set(findobj('tag','PTSpanel'),'title','Point type selection',...
        'fontname',font,'fontsize',8);   
    set(findobj('tag','GMaRSSpanel'),'title',...
        'Geopotential model and reference system selection','fontname',...
        font,'fontsize',8); 
    set(findobj('tag','CPaOSpanel'),'title',...
        'Calculated parameters and output selection','fontname',...
        font,'fontsize',8);

else  
        
    switch(vstpar)
        
        case('import_GGM') %Click on the Browse... button in the Geopotential model and reference system selection panel
            
                [GGMname,GGMadresar]=uigetfile('*.*','Select GGM File');
                if GGMname==0
                else
                    set(findobj('tag','R'),'userdata',GGMname);
                    set(findobj('tag','ell'),'userdata',GGMadresar);
                    
                    set(findobj('tag','nameGGM'),'string',GGMname); %Display the name of the imported GGM file
                end
                
        case('input') %Click on the Browse... button in the Point type selection panel
            
            [loadname,loadadresar]=uigetfile('*.*','Select File Containing Computational Points');
            if loadname==0
            else
                set(findobj('tag','use'),'userdata',loadname);
                set(findobj('tag','diskcheck'),'userdata',loadadresar);
                
                set(findobj('tag','LoadFile'),'string',loadname); %Display the name of the imported GGM file
            end
        
        case('Output_folder') %Click on the Output folder and file button
            
            [outname,outadresar]=uiputfile('*.*');
            if outname==0
            else
                if find(outname=='.')>0
                    outname=outname(1:(find(outname=='.', 1,'last')-1));
                end
                set(findobj('tag','R_text'),'userdata',outname);
                set(findobj('tag','ell_text'),'userdata',outadresar);
                
                set(findobj('tag','OutFile'),'string',outname); %Display name of the output file
            end
                        
        case('Display_data_settings') %Click on the Display data settings
            
            %Main window
            D=figure('units','pixels','numbertitle','off','name',...
                'Display data settings','color',[R1 G1 B1],'position',...
                [screenres(3)/2-550/2 screenres(4)/2-450/2 550 450],'tag','oknoDDS','menubar','none');

            %Panel
            uipanel('Units','normalized','position',[0 0 0.28 1],...
                'backgroundcolor',[R1 G1 B1],'HighlightColor',[R2 G2 B2],...
                'fontsize',8,'fontname',font);
            
            %Display data checkbox
            uicontrol('Units','normalized','position',[0.02 0.92 0.2 0.05],...
                'style','checkbox','string','Display data','backgroundcolor',...
                [R1 G1 B1],'tag','Display');
            
            display_data=get(findobj('tag','DDS'),'userdata');
            if display_data==0
            elseif display_data==1
                set(findobj('tag','Display'),'value',1);
            end

            %Text next to the Display data checkbox
            g0=uicontrol('Units','normalized','position',[0.31 0.85 0.67 0.13],...
                'style','text','backgroundcolor',[R1 G1 B1],'string',...
                'To plot the results, select this checkbox. This option works only if "Point type selection" is set to "Grid".');
            set(g0,'HorizontalAlignment','left');
            
            %Mapping toolbox vs. basic Matlab function imagesc checkbox
            uicontrol('Units','normalized','position',[0.02 0.79 0.25 0.05],...
                'style','checkbox','string','Mapping Toolbox','backgroundcolor',...
                [R1 G1 B1],'tag','Mapp_tool');
            
            Display_data=get(findobj('tag','h_string'),'userdata');
            if Display_data==0
            elseif Display_data==1
                set(findobj('tag','Mapp_tool'),'value',1);
            end
            
            %Text next to the Display data checkbox
            g01=uicontrol('Units','normalized','position',[0.31 0.72 0.67 0.16],...
                'style','text','backgroundcolor',[R1 G1 B1],'string',...
                'Select this checkbox if Mapping Toolbox is to be used to plot the data (may be painfully slow). Leave this checkbox unselected if the imagesc function is to be used (does not require Mapping Toolbox and is substantially faster).');
            set(g01,'HorizontalAlignment','left');
            
            %Text Graphic format
            uicontrol('Units','normalized','position',[0.04 0.65 0.2 0.05],...
                'style','text','string','Graphic format','backgroundcolor',...
                [R1 G1 B1],'tag','format');
            
            %Graphic format pop-up menu
            uicontrol('Units','normalized','position',...
                [0.06 0.55 0.16 0.1],'style','popup','string',...
                '*.bmp|*.emf|*.eps|*.jpeg|*.pdf|*.png|*.tiff',...
                'backgroundcolor',[R2 G2 B2],'tag','pripona');
            
            prip=get(findobj('tag','nmin'),'userdata');     
            if isempty(prip)
                set(findobj('tag','pripona'),'value',6);
            else  
                set(findobj('tag','pripona'),'value',prip);
            end
            
            %Text next to the Graphic format file
            g1=uicontrol('Units','normalized','position',[0.31 0.57 0.66 0.13],...
                'style','text','backgroundcolor',[R1 G1 B1],'string',...
                'Select one of the graphic format files. For a vector output it is recommended to use *.eps graphic file and *.png format for a bitmap output.');
            set(g1,'HorizontalAlignment','left');
            
            %Text Colormap
            uicontrol('Units','normalized','position',[0.06 0.505 0.16 0.05],...
                'style','text','string','Colormap',...
                'backgroundcolor',[R1 G1 B1]);
            
            %Colormap pop-up menu
            uicontrol('Units','normalized','position',[0.06 0.405 0.16 0.1],...
                'style','popup','string',...
                'jet|HSV|hot|cool|spring|summer|autumn|winter|gray|bone|copper|pink|lines',...
                'backgroundcolor',[R2 G2 B2],'tag','colormap');
            set(findobj('tag','colormap'),'value',2);
            
            color=get(findobj('tag','nmax'),'userdata');     
            if isempty(color)
                set(findobj('tag','colormap'),'value',1);
            else  
                set(findobj('tag','colormap'),'value',color);
            end
            
            %Text next to the colormap
            g1=uicontrol('Units','normalized','position',[0.31 0.355 0.65 0.2],...
                'style','text','backgroundcolor',[R1 G1 B1],'string',...
                'Select a colormap of the output file. Frequently, the jet colormap is used, which ranges from blue to red, and passes through the colors cyan, yellow, and orange.');
            set(g1,'HorizontalAlignment','left');
            
            %Text Number of colors
            uicontrol('Units','normalized','position',[0.04 0.34 0.21 0.05],...
                'style','text','string','Number of colors',...
                'backgroundcolor',[R1 G1 B1]);
            
            %Value of number of colors
            uicontrol('Units','normalized','position',[0.09 0.28 0.1 0.06],...
                'style','edit','string',15,'backgroundcolor',...
                [R2 G2 B2],'tag','skala');
            
            %Text next to the number of colors
            g1=uicontrol('Units','normalized','position',[0.31 0.255 0.65 0.15],...
                'style','text','backgroundcolor',[R1 G1 B1]);
            set(g1,'HorizontalAlignment','left','string',...
                'Enter a number of colors of the selected colormap. If the the Mapping Toolbox checkbox is selected and a larger number of colors is entered, the processing time may significantly increase.');
            
            ncolor=get(findobj('tag','text_nmin'),'userdata');
            if isempty(ncolor)                
            else
                set(findobj('tag','skala'),'string',ncolor);
            end
            
            %Text DPI 
            uicontrol('Units','normalized','position',[0.09 0.20 0.1 0.05],...
                'style','text','string','DPI','backgroundcolor',[R1 G1 B1]);
            
            %Value of DPI
            uicontrol('Units','normalized','position',[0.09 0.15 0.1 0.06],...
                'style','edit','string',300,'backgroundcolor',[R2 G2 B2],...
                'tag','DPI');
            
            DPI=get(findobj('tag','text_nmax'),'userdata');
            if isempty(DPI)                
            else
                set(findobj('tag','DPI'),'string',DPI);
            end
            
            %Text next to the DPI
            g1=uicontrol('Units','normalized','position',[0.31 0.00 0.65 0.2],...
                'style','text','backgroundcolor',[R1 G1 B1],'string',...
                'Enter a value of dots per inch of the output file.');
            set(g1,'HorizontalAlignment','left');
                        
            %OK button
            uicontrol('Units','normalized','position',[0.35 0.02 0.13 0.08],...
                'style','pushbutton','string','OK','Callback','GrafLab OKDDS');
    
            %Close button
            uicontrol('Units','normalized','position',[0.55 0.02 0.13 0.08],...
                'style','pushbutton','string','Close','Callback',...
                'GrafLab CloseDDS');
            
            %setting font to Cambria
            set(get(D,'children'),'fontname',font,'fontsize',10)
            
        case('fnALFs') %Click on the Computation of fnALFs
            
            %Main window
            F=figure('units','pixels','numbertitle','off','name',...
                'Computation of fully normalized associated Legendre functions','color',...
                [R1 G1 B1],'position',[screenres(3)/2-550/2 screenres(4)/2-350/2 550 350],'tag',...
                'oknofnALFs','menubar','none'); 
            
            %Radio button group
            u0 = uibuttongroup('visible','on','units','normalized',...
                'Position',[0 0 .26 1],'tag','volbaALFs'); 

            %The first radiobutton
            u1 = uicontrol('units','normalized','Style','Radio','pos',...
                [0.1 0.7 0.9 0.3],'parent',u0,'HandleVisibility','on',...
                'tag','rbutton1'); 
            set(u1,'String','<html>Standard forward<br>column method',...
                'fontname',font,'fontsize',10);
            note1=uicontrol('Units','normalized','position',...
                [0.3 0.72 0.66 0.24],'style','text','backgroundcolor',...
                [R1 G1 B1]);
            set(note1,'HorizontalAlignment','left','string',...
                'It is recommended to use the standard forward column method for all latitudes up to the maximum degree 1800. However, this method may also be used for the latitudes <0 deg, 56 deg> and <80 deg, 90 deg> up to the maximum degree 2190.');
            
            %The second radiobutton
            u2 = uicontrol('units','normalized','Style','Radio','pos',...
                [0.1 0.39 0.9 0.3],'parent',u0,'HandleVisibility','on',...
                'tag','rbutton2'); 
            set(u2,'String',...
                '<html>Modified forward<br>column method<br>combined with<br>Horner''s scheme',...
                'fontname',font,'fontsize',10);
            note2=uicontrol('Units','normalized','position',[0.3 0.38 0.65 0.3],...
                'style','text','backgroundcolor',[R1 G1 B1]);
            set(note2,'HorizontalAlignment','left','string',...
                'It is recommended to use the modified forward column method combined with Horner''s scheme for all latitudes and maximum degrees ranging from 1801 to 2700. This method may also be used for lower degrees than 1801, but cannot be applied to degrees higher than 2700 due to the overflow problem.');            
            
            %The third radiobutton
            u3 = uicontrol('units','normalized','Style','Radio','pos',...
                [0.1 0.12 0.9 0.3],'parent',u0,'HandleVisibility','on',...
                'tag','rbutton3'); 
            set(u3,'String','<html>Extended-range<br>arithmetic',...
                'fontname',font,'fontsize',10);
            note3=uicontrol('Units','normalized','position',...
                [0.3 0.018 0.65 0.3],'style','text','backgroundcolor',...
                [R1 G1 B1]);
            set(note3,'HorizontalAlignment','left','string',...
                'The extended-range arithmetic approach may be used for all latitudes up to an arbitrary degree essentially.');                        
            
            %OK button
            uicontrol('Units','normalized','position',[0.35 0.02 0.13 0.08],...
                'style','pushbutton','string','OK','Callback','GrafLab OKfnALFs');
    
            %Close button
            uicontrol('Units','normalized','position',[0.55 0.02 0.13 0.08],...
                'style','pushbutton','string','Close','Callback','GrafLab ClosefnALFs');
            
            %The chosen approach for computation of fnALFs
            volbaALFs=get(findobj('tag','fnALFs'),'userdata');
            if isempty(volbaALFs)
                volbaALFs=1;
            end
            
            %Mark the chosen radiobutton
            if volbaALFs==1
                set(u1,'value',1);
            elseif volbaALFs==2
                set(u2,'value',1);
            elseif volbaALFs==3
                set(u3,'value',1);
            end
            
            set(get(F,'children'),'fontname',font,'fontsize',10);
            
        case('ClosefnALFs') %Click on the Close button in the Computation of fnALFs window
            
            close
            
        case('OKfnALFs') %Click on the OK button in the Computation of fnALFs window
            
            volbaALFs=get(findobj('tag','volbaALFs'),'selectedobject');

            if get(findobj('tag','rbutton1'),'value')==1
                volbaALFs=1;
            elseif get(findobj('tag','rbutton2'),'value')==1
                volbaALFs=2;
            elseif get(findobj('tag','rbutton3'),'value')==1
                volbaALFs=3;
            end
            
            set(findobj('tag','fnALFs'),'userdata',volbaALFs);            
            
            close
            
        case('OKDDS') %Click on the OK button in the Display data settings window
            
            display_data=get(findobj('tag','Display'),'value');
            set(findobj('tag','DDS'),'userdata',display_data);   
            
            Display_data=get(findobj('tag','Mapp_tool'),'value');
            set(findobj('tag','h_string'),'userdata',Display_data);
            
            prip=get(findobj('tag','pripona'),'value');
            set(findobj('tag','nmin'),'userdata',prip);
            
            color=get(findobj('tag','colormap'),'value');
            set(findobj('tag','nmax'),'userdata',color);

            ncolor=get(findobj('tag','skala'),'string');
            ncolor=str2double(ncolor);
            set(findobj('tag','text_nmin'),'userdata',ncolor);
            
            DPI=get(findobj('tag','DPI'),'string');
            DPI=str2double(DPI);
            set(findobj('tag','text_nmax'),'userdata',DPI);
            
            %Check of the entered value of DPI and number of colors
            if isnan(DPI)==1 
                errordlg('The DPI value must be a number.',...
                    'Display data settings');
                error('The DPI value must be a number.');
            end
            
            if isnan(ncolor)==1
                errordlg('The entered value of number of colors is not correct.',...
                    'Display data settings');
                error('The entered value of number of colors is not correct.');
            end
            
            if rem(DPI,1)~=0 || DPI<0
                errordlg('The DPI value must be a positive integer.',...
                    'Display data settings');
                error('The DPI value must be a positive integer.');
            end
                   
            if rem(ncolor,1)~=0 || ncolor<2
                errordlg('The value of number of colors must be an integer larger than 1.',...
                    'Display data settings');
                error('The value of number of colors must be an integer larger than 1.');
            end
            
            close
            
        case('CloseDDS') %Click on the Close button in the Display data settings window
            
            close
            
        case('OK') %Click on the OK button in the main GrafLab window 
        
            set(findobj('tag','hlasky'),'string',...
                    '','fontsize',8,'foregroundcolor','k'); drawnow;   
            
            if nargin==1 %Working with the GUI
                GUI=1;
                Status_bar=0;
            else
                GUI=0; %Working without the GUI
            end
                
            if GUI==1 %If working with the GUI
                
                %Identification of point type selection
                volbagridcheck=get(findobj('tag','rbutton1type'),'value');
                volbaloadcheck=get(findobj('tag','rbutton2type'),'value');
                volbadiskcheck=get(findobj('tag','rbutton3type'),'value');

                %GM and R of GGM
                GM=str2double(get(findobj('tag','GM'),'string'));            
                R=str2double(get(findobj('tag','R'),'string'));
                
                %Reference ellipsoid
                ellipsoid=get(findobj('tag','ell'),'value');
                
                %Identification of type of the input coordinates
                coord=get(findobj('tag','rbutton2coord'),'value');
                
                %Identification of computation of fnALFs approach
                volbaALFs=get(findobj('tag','fnALFs'),'userdata');
                
                %nmin value
                nmin=str2double(get(findobj('tag','nmin'),'string'));
                
                %nmax value
                nmax=str2double(get(findobj('tag','nmax'),'string'));
                
                %Use maximum degree of the GGM
                use_nmax_GGM=get(findobj('tag','use'),'value');
                
                %Identification of Display data
                display_data=get(findobj('tag','DDS'),'userdata');
                if display_data==1
                    Display_data=get(findobj('tag','h_string'),'userdata');
                    if Display_data==0 
                        display_data=2;
                    elseif Display_data==1
                        display_data=1;
                    end
                    clear Display_data
                end
                ncolor=get(findobj('tag','text_nmin'),'userdata');
                DPI=get(findobj('tag','text_nmax'),'userdata');
                volbaformat=get(findobj('tag','nmin'),'userdata');
                colmap=get(findobj('tag','nmax'),'userdata');
                
                %Identification of computation of commission error
                STD=get(findobj('tag','STD'),'value');
                
                %Import of data file containing ellipsoidal coordinates
                loadname=get(findobj('tag','use'),'userdata');
                loadadresar=get(findobj('tag','diskcheck'),'userdata');
                
                %Selection of functionals of the geopotential
                volbapar1=get(findobj('tag','P1'),'value'); %ID of the first functional
                volbapar2=get(findobj('tag','P2'),'value'); %ID of the second functional
                volbapar3=get(findobj('tag','P3'),'value'); %ID of the third functional
                volbapar4=get(findobj('tag','P4'),'value'); %ID of the fourth functional
                volbapar=[volbapar1;volbapar2;volbapar3;volbapar4];
                
                %Identification of the output folder and file
                outname=get(findobj('tag','R_text'),'userdata');
                outadresar=get(findobj('tag','ell_text'),'userdata');
                
                %Identification of the input GGM file
                GGMname=get(findobj('tag','R'),'userdata');
                GGMadresar=get(findobj('tag','ell'),'userdata');
                
                %Export report
                Export_report=get(findobj('tag','report'),'value');
                
                %Export data in "mat" format
                Export_data_mat=get(findobj('tag','datamat'),'value');
                
                %Export data in "txt" format
                Export_data_txt=get(findobj('tag','export'),'value');
                
                %Import of the variance-covariance matrix file                           
                GGMcovname=get(findobj('tag','R'),'userdata');
                GGMcovadresar=get(findobj('tag','ell'),'userdata');
                
            elseif GUI==0 %If working without the GUI
            
                volbapar=Functional(:);
                clear Functional
                if length(volbapar)==1 %1 functional to be computed
                    volbapar1=volbapar(1); %ID of the first functional
                    volbapar2=1; volbapar3=1; volbapar4=1;
                elseif length(volbapar)==2 %2 functionals to be computed
                    volbapar1=volbapar(1); %ID of the first functional
                    volbapar2=volbapar(2); %ID of the second functional
                    volbapar3=1; volbapar4=1;
                elseif length(volbapar)==3 %3 functionals to be computed
                    volbapar1=volbapar(1); %ID of the first functional
                    volbapar2=volbapar(2); %ID of the second functional
                    volbapar3=volbapar(3); %ID of the third functional
                    volbapar4=1;
                elseif length(volbapar)==4 %4 functionals to be computed
                    volbapar1=volbapar(1); %ID of the first functional
                    volbapar2=volbapar(2); %ID of the second functional
                    volbapar3=volbapar(3); %ID of the third functional
                    volbapar4=volbapar(4); %ID of the fourth functional
                else
                    error('Please choose one to four functionals to compute.')
                end
                
                %Some error checks of the "point_type" variable
                if isempty(point_type) || length(point_type)>1 || ~isnumeric(point_type) || rem(point_type,1)~=0
                    error('The variable "point_type" must be an integer (0 -- Grid, 1 -- Load data, 2 -- Point-wise).')
                end
            
                if point_type==0 || point_type==1 || point_type==2
                else
                    error('The variable "point_type" must be an integer (0 -- Grid, 1 -- Load data, 2 -- Point-wise).')
                end
                
                %Point type selection
                if point_type==0 %Grid-wise computation
                    volbagridcheck=1; %Grid
                    volbaloadcheck=0; %Load
                    volbadiskcheck=0; %Point-wise
                elseif point_type==1 %Load data
                    volbagridcheck=0; %Grid
                    volbaloadcheck=1; %Load
                    volbadiskcheck=0; %Point-wise
                elseif point_type==2 %Point-wise computation
                    volbagridcheck=0; %Grid
                    volbaloadcheck=0; %Load
                    volbadiskcheck=1; %Point-wise
                end
                
                %Some error checks of the "Functional_or_Commission" variable
                if isempty(Functional_or_Commission) || length(Functional_or_Commission)>1 || ~isnumeric(Functional_or_Commission) || rem(point_type,1)~=0
                    error('The variable "Functional_or_Commission" must be an integer (0 -- functional, 1 -- commission error).')
                end
            
                if Functional_or_Commission==0 || Functional_or_Commission==1
                else
                    error('The variable "Functional_or_Commission" must be an integer (0 -- functional, 1 -- commission error).')
                end
                
                %Computation of functionals or commission errors
                STD=Functional_or_Commission; clear Functional_or_Commission
                
                %Approach to the computation of Legendre functions
                volbaALFs=fnALFs; clear fnALFs
                
                if isempty(Display_data)
                    Display_data=0;
                end
                
                %Some error checks of the "Display_data" variable
                if length(Display_data)>1 || ~isnumeric(Display_data)
                    error('The variable "Display_data" must be an integer (0 -- do not display data, 1 -- display data using Mapping Toolbox, 2 -- display data using the "imagesc" function.')
                end
            
                if Display_data==2 || Display_data==1 || Display_data==0
                else
                    error('The variable "Display_data" must be an integer (0 -- do not display data, 1 -- display data using Mapping Toolbox, 2 -- display data using the "imagesc" function).')
                end
                
                if Display_data==1 || Display_data==2
                    if isempty(Graphic_format) || any(rem(Graphic_format,1)~=0) || any(Graphic_format<1) || any(Graphic_format>7)
                        error('When displaying data, the variable "Graphic_format" must be an integer: 1, ..., 7 (1 - bmp, 2 - emf, 3 - eps, 4 - jpeg, 5 - pdf, 6 - png, 7 - tiff).')
                    end
                    
                    if isempty(Colormap) || any(rem(Colormap,1)~=0) || any(Colormap<1) || any(Colormap>13)
                        error('When displaying data, the variable "Colormap" must be an integer: 1, ..., 13 (1 - jet, 2 - hsv, 3 - hot, 4 - cool, 5 - spring, 6 - summer, 7 - autumn, 8 - winter, 9 - gray, 10 - bone, 11 - copper, 12 - pink, 13 - lines).')
                    end
                    
                    if isempty(Number_of_colors) || length(Number_of_colors)>1 || ~isnumeric(Number_of_colors) || rem(Number_of_colors,1)~=0 || Number_of_colors<=0
                        error('When displaying data, the variable "Number_of_colors" must be a positive integer.')
                    end
                    
                    if isempty(DPI) || length(DPI)>1 || ~isnumeric(DPI) || rem(DPI,1)~=0 || DPI<=0
                        error('When displaying data, the variable "DPI" must be a positive integer.')
                    end
                end
                
                %Display data settings
                display_data=Display_data; clear Display_data
                ncolor=Number_of_colors; clear Number_of_colors
                colmap=Colormap; clear Colormap
                volbaformat=Graphic_format; clear Graphic_format
                
                %Identification of the name and the path to the output file
                if isempty(find(Output_path=='/' | Output_path=='\', 1, 'last' ))
                    outname=Output_path;
                    outadresar=[];
                else
                    temp=find(Output_path=='/' | Output_path=='\');
                    outname=Output_path(max(temp)+1:end);
                    outadresar=Output_path(1:max(temp));
                    clear temp
                end
                
                %Identification of the name and the path to the input GGM
                %file
                if STD==0
                    if isempty(find(GGM_path=='/' | GGM_path=='\', 1, 'last' ))
                        GGMname=GGM_path;
                        GGMadresar=[];
                    else
                        temp=find(GGM_path=='/' | GGM_path=='\');
                        GGMname=GGM_path(max(temp)+1:end);
                        GGMadresar=GGM_path(1:max(temp));
                        clear temp
                    end
                elseif STD==1
                    if isempty(find(GGM_path=='/' | GGM_path=='\', 1, 'last' ))
                        GGMcovname=GGM_path;
                        GGMcovadresar=[];
                    else
                        temp=find(GGM_path=='/' | GGM_path=='\');
                        GGMcovname=GGM_path(max(temp)+1:end);
                        GGMcovadresar=GGM_path(1:max(temp));
                        clear temp
                    end
                end
                
                %Identification of the name and the path to the DTM file
                if isempty(find(DTM_path=='/' | DTM_path=='\', 1, 'last' ))
                    loadnameDMR=DTM_path;
                    loadadresarDMR=[];
                else
                    temp=find(DTM_path=='/' | DTM_path=='\');
                    loadnameDMR=DTM_path(max(temp)+1:end);
                    loadadresarDMR=DTM_path(1:max(temp));
                    clear temp
                end
                
                %Identification of the name and the path to the input data
                %file (load data option)
                if isempty(find(Input_data_path=='/' | Input_data_path=='\', 1, 'last' ))
                    loadname=Input_data_path;
                    loadadresar=[];
                else
                    temp=find(Input_data_path=='/' | Input_data_path=='\');
                    loadname=Input_data_path(max(temp)+1:end);
                    loadadresar=Input_data_path(1:max(temp));
                    clear temp
                end
                
                %Names of the functionals (to be displayed in the report file,
                %in names of the displayed functionals, etc.)
                P1={'' 'Deflection_of_the_vertical_eta' 'Deflection_of_the_vertical_xi' 'Deflection_of_the_vertical_Theta' 'Disturbing_potential' 'Disturbing_tensor_Trr_Tpp_Tll' 'Disturbing_tensor_Trp_Trl_Tpl' 'Disturbing_tensor_Txx_Tyy_Tzz' 'Disturbing_tensor_Txy_Txz_Tyz' 'Geoid_undulation' 'Gravitational_potential' 'Gravitational_tensor_Vrr_Vpp_Vll' 'Gravitational_tensor_Vrp_Vrl_Vpl' 'Gravitational_tensor_Vxx_Vyy_Vzz' 'Gravitational_tensor_Vxy_Vxz_Vyz' 'Gravity_vector_gX_gY_gZ' 'Gravity_sa' 'Gravity_potential' 'Gravity_anomaly_sa' 'Gravity_disturbance' 'Gravity_disturbance_sa' 'Height_anomaly_ell' 'Height_anomaly' '2nd_radial_der_of_disturbing_potential' '2nd_radial_der_of_gravity_potential'};
                P2={'' 'Deflection_of_the_vertical_eta' 'Deflection_of_the_vertical_xi' 'Deflection_of_the_vertical_Theta' 'Disturbing_potential' 'Disturbing_tensor_Trr_Tpp_Tll' 'Disturbing_tensor_Trp_Trl_Tpl' 'Disturbing_tensor_Txx_Tyy_Tzz' 'Disturbing_tensor_Txy_Txz_Tyz' 'Geoid_undulation' 'Gravitational_potential' 'Gravitational_tensor_Vrr_Vpp_Vll' 'Gravitational_tensor_Vrp_Vrl_Vpl' 'Gravitational_tensor_Vxx_Vyy_Vzz' 'Gravitational_tensor_Vxy_Vxz_Vyz' 'Gravity_vector_gX_gY_gZ' 'Gravity_sa' 'Gravity_potential' 'Gravity_anomaly_sa' 'Gravity_disturbance' 'Gravity_disturbance_sa' 'Height_anomaly_ell' 'Height_anomaly' '2nd_radial_der_of_disturbing_potential' '2nd_radial_der_of_gravity_potential'};
                P3={'' 'Deflection_of_the_vertical_eta' 'Deflection_of_the_vertical_xi' 'Deflection_of_the_vertical_Theta' 'Disturbing_potential' 'Disturbing_tensor_Trr_Tpp_Tll' 'Disturbing_tensor_Trp_Trl_Tpl' 'Disturbing_tensor_Txx_Tyy_Tzz' 'Disturbing_tensor_Txy_Txz_Tyz' 'Geoid_undulation' 'Gravitational_potential' 'Gravitational_tensor_Vrr_Vpp_Vll' 'Gravitational_tensor_Vrp_Vrl_Vpl' 'Gravitational_tensor_Vxx_Vyy_Vzz' 'Gravitational_tensor_Vxy_Vxz_Vyz' 'Gravity_vector_gX_gY_gZ' 'Gravity_sa' 'Gravity_potential' 'Gravity_anomaly_sa' 'Gravity_disturbance' 'Gravity_disturbance_sa' 'Height_anomaly_ell' 'Height_anomaly' '2nd_radial_der_of_disturbing_potential' '2nd_radial_der_of_gravity_potential'};
                P4={'' 'Deflection_of_the_vertical_eta' 'Deflection_of_the_vertical_xi' 'Deflection_of_the_vertical_Theta' 'Disturbing_potential' 'Disturbing_tensor_Trr_Tpp_Tll' 'Disturbing_tensor_Trp_Trl_Tpl' 'Disturbing_tensor_Txx_Tyy_Tzz' 'Disturbing_tensor_Txy_Txz_Tyz' 'Geoid_undulation' 'Gravitational_potential' 'Gravitational_tensor_Vrr_Vpp_Vll' 'Gravitational_tensor_Vrp_Vrl_Vpl' 'Gravitational_tensor_Vxx_Vyy_Vzz' 'Gravitational_tensor_Vxy_Vxz_Vyz' 'Gravity_vector_gX_gY_gZ' 'Gravity_sa' 'Gravity_potential' 'Gravity_anomaly_sa' 'Gravity_disturbance' 'Gravity_disturbance_sa' 'Height_anomaly_ell' 'Height_anomaly' '2nd_radial_der_of_disturbing_potential' '2nd_radial_der_of_gravity_potential'};
                
                if isempty(Status_bar)
                    Status_bar=0;
                end
                
                if Status_bar==0 || Status_bar==1
                else
                    error('The variable "Status_bar" can take the value 0 (off), 1 (on) or can be left empty (off).')
                end
                
            end
            
            %Error messages for point type selection
            if volbagridcheck==0 && volbadiskcheck==0 && volbaloadcheck==0
                if GUI==1 %If working with the GUI
                    errordlg('Please select a point type.','Point type selection error');
                end
                error('Please select a point type.');
            elseif volbagridcheck==1 && volbadiskcheck==1 && volbaloadcheck==1
                if GUI==1 %If working with the GUI
                    errordlg('Please select one point type only.',...
                        'Error in point type selection');
                end
                error('Please select one point type only.');
            elseif volbagridcheck==1 && volbadiskcheck==1
                if GUI==1 %If working with the GUI
                    errordlg('Please select one point type only.',...
                        'Error in point type selection');
                end
                error('Please select one point type only.');
            elseif volbadiskcheck==1 && volbaloadcheck==1
                if GUI==1 %If working with the GUI
                    errordlg('Please select one point type only.',...
                        'Error in point type selection');
                end
                error('Please select one point type only.');
            elseif volbagridcheck==1 && volbaloadcheck==1
                if GUI==1 %If working with the GUI
                    errordlg('Please select one point type only.',...
                        'Error in point type selection');
                end
                error('Please select one point type only.');
            end            
                            
            %Ellipsoid
            if length(ellipsoid)==1
                if ellipsoid==1 %GRS80              
                    GMEl=3986005*10^8; %Geocentric gravitational constant of GRS80
                    aEl=6378137; %Semimajor axis of GRS80
                    eEl=sqrt(0.006694380022903416); %First eccentricity of GRS80
                    omegaEl=7292115*10^-11; %Angular velocity of GRS80
                    CEl_20=-108263*10^-8/sqrt(5); %Fully normalized C_20 of GRS80
                elseif ellipsoid==2 %WGS84
                    GMEl=3986004.418*10^8; %Geocentric gravitational constant of WGS84
                    aEl=6378137; %Semimajor axis of WGS84
                    fEl=1/298.257223563; %Flattening of WGS84
                    omegaEl=7292115*10^-11; %Angular velocity of WGS84
                    CEl_20=-0.484166774985*10^-3; %Fully normalized C_20 of WGS84
                    eEl=sqrt(fEl*(2-fEl)); %First eccentricity of WGS84
                else
                    error('The variable "ellipsoid" must be either a scalar (1 -- GRS80, 2 -- WGS84) or a vector with 5 elements [GMEl aEl eEl CEl_20 omegaEl].')
                end
            elseif length(ellipsoid)==5
                GMEl=ellipsoid(1); %Geocentric gravitational constant
                aEl=ellipsoid(2); %Semimajor axis
                eEl=ellipsoid(3); %First eccentricity
                CEl_20=ellipsoid(4); %Fully normalized C_20
                omegaEl=ellipsoid(5); %Angular velocity
            else
                error('The variable "ellipsoid" must be either a scalar (1 -- GRS80, 2 -- WGS84) or a vector with 5 elements [GMEl aEl eEl CEl_20 omegaEl].')
            end
            
            %Some error checks of the GM value
            if isempty(GM) || ~isnumeric(GM) || length(GM)>1
                if GUI==1
                    errordlg('The value of "GM" must a real-valued scalar.',...
                        'Error in geopotential model and reference system selection')
                end
                error('The value of "GM" must a real-valued scalar.')
            end
            if GM<=0
                if GUI==1 %Working with the GUI
                    errordlg('The value of "GM" must be larger than zero.',...
                        'Error in geopotential model and reference system selection')
                end
                error('The value of "GM" must be larger than zero.')
            end
            if isnan(GM)
                if GUI==1 %Working with the GUI
                    errordlg('Please input the value of "GM".',...
                        'Error in geopotential model and reference system selection')
                end
                error('Please input the value of "GM".')
            end
            
            %Some error checks of the R value
            if isempty(R) || ~isnumeric(R) || length(R)>1
                if GUI==1
                    errordlg('The value of "R" must a real-valued scalar.',...
                        'Error in geopotential model and reference system selection')
                end
                error('The value of "R" must a real-valued scalar.')
            end
            if R<=0
                if GUI==1 %Working with the GUI
                    errordlg('The value of "R" must be larger than zero.',...
                        'Error in geopotential model and reference system selection')
                end
                error('The value of "R" value must be larger than zero.')
            end     
            if isnan(R)
                if GUI==1 %Working with the GUI
                    errordlg('Please input the value of "R".',...
                        'Error in geopotential model and reference system selection')
                end
                error('Please input the value of "R".')
            end        
            
            %Some error checks of the "coord" variable
            if isempty(coord) || length(coord)>1 || ~isnumeric(coord) || rem(coord,1)~=0
                error('The variable "coord" must be an integer (0 -- ellipsoidal coordinates, 1 -- spherical coordinates).')
            end
            
            if coord==1 || coord==0
            else
                error('The variable "coord" must be an integer (0 -- ellipsoidal coordinates, 1 -- spherical coordinates).')
            end
            
            %If the approach to compute fnALFs is not specified, GrafLab
            %automatically uses the standard forward column method
            if isempty(volbaALFs)
                volbaALFs=1;
            end
            
            %Some error checks of the "volbaALFs" ("fnALFs") variable
            if isempty(volbaALFs) || length(volbaALFs)>1 || ~isnumeric(volbaALFs) || rem(volbaALFs,1)~=0
                error('The variable "fnALFs" must be an integer (1 - standard forward column method, 2 - modified forward column method, 3 - extended-range arithmetic).')
            end
            
            if volbaALFs==1 || volbaALFs==2 || volbaALFs==3
            else
                error('The variable "fnALFs" must be an integer (1 - standard forward column method, 2 - modified forward column method, 3 - extended-range arithmetic).')
            end
            
            %Total number of functionals to be computed
            pocetpar=length(find(volbapar>1));

            %Error messages for selection of functionals of the geopotential
            if pocetpar==1
                if volbapar1>1 && volbapar2==1 && volbapar3==1 && volbapar4==1
                else
                    if GUI==1 %If working with the GUI
                        errordlg('Please select the functional of the geopotential in the first pop-up menu.',...
                            'Calculated parameters and output selection');
                    end
                    error('Please select the functional of the geopotential in the first pop-up menu.');
                end
            elseif pocetpar==2
                if volbapar1>1 && volbapar2>1 && volbapar3==1 && volbapar4==1
                    if volbapar1 == volbapar2
                        if GUI==1 %If working with the GUI
                            errordlg('An identical functional of the geopotential cannot be computed simultaneously more than once.',...
                                'Calculated parameters and output selection');
                        end
                        error('An identical functional of the geopotential cannot be computed simultaneously more than once.');
                    end
                else
                    if GUI==1 %If working with the GUI
                        errordlg('Please select the functionals of the geopotential in the first and the second pop-up menu.',...
                            'Calculated parameters and output selection');
                    end
                    error('Please select the functionals of the geopotential in the first and the second pop-up menu.');
                end
                
                if length(nonzeros(volbapar==10 | volbapar==23))==1 %If geoid or height anomaly is to be computed
                    if GUI==1 %If working with the GUI
                        errordlg('Geoid_undulation and/or Height_anomaly cannot be computed simultaneously with other functionals. However, Geoid_undulation and Height_anomaly alone can be computed simultaneously.',...
                            'Calculated parameters and output selection');
                    end
                    error('Geoid_undulation and/or Height_anomaly cannot be computed simultaneously with other functionals. However, Geoid_undulation and Height_anomaly alone can be computed simultaneously.');
                end
            elseif pocetpar==3
                if volbapar1>1 && volbapar2>1 && volbapar3>1 && volbapar4==1
                    if (volbapar1 == volbapar2) || (volbapar1 == volbapar3) || (volbapar2 == volbapar3)
                        if GUI==1 %If working with the GUI
                            errordlg('An identical functional of the geopotential cannot be computed simultaneously more than once.',...
                                'Calculated parameters and output selection');
                        end
                        error('An identical functional of the geopotential cannot be computed simultaneously more than once.');
                    end
                else
                    if GUI==1 %If working with the GUI
                        errordlg('Please select the functionals of the geopotential in the first, second and third pop-up menu.',...
                            'Calculated parameters and output selection');
                    end
                    error('Please select the functionals of the geopotential in the first, second and third pop-up menu.');
                end
                
                if length(nonzeros(volbapar==10 | volbapar==23))==1 || length(nonzeros(volbapar==10 | volbapar==23))==2 %If geoid or height anomaly is to be computed
                    if GUI==1 %If working with the GUI
                        errordlg('Geoid_undulation and/or Height_anomaly cannot be computed simultaneously with other functionals. However, Geoid_undulation and Height_anomaly alone can be computed simultaneously.',...
                            'Calculated parameters and output selection');
                    end
                    error('Geoid_undulation and/or Height_anomaly cannot be computed simultaneously with other functionals. However, Geoid_undulation and Height_anomaly alone can be computed simultaneously.');
                end
            elseif pocetpar==4
                if (volbapar1 == volbapar2) || (volbapar1 == volbapar3) || (volbapar1 == volbapar4) || (volbapar2 == volbapar3) || (volbapar2 == volbapar4) || (volbapar3 == volbapar4)
                    if GUI==1 %If working with the GUI
                        errordlg('An identical functional of the geopotential cannot be computed simultaneously more than once.',...
                            'Calculated parameters and output selection');
                    end
                    error('An identical functional of the geopotential cannot be computed simultaneously more than once.');
                end
                
                if length(nonzeros(volbapar==10 | volbapar==23))==1 || length(nonzeros(volbapar==10 | volbapar==23))==2 || length(nonzeros(volbapar==10 | volbapar==23))==3 %If geoid or height anomaly is to be computed
                    if GUI==1 %If working with the GUI
                        errordlg('Geoid_undulation and/or Height_anomaly cannot be computed simultaneously with other functionals. However, Geoid_undulation and Height_anomaly alone can be computed simultaneously.',...
                            'Calculated parameters and output selection');
                    end
                    error('Geoid_undulation and/or Height_anomaly cannot be computed simultaneously with other functionals. However, Geoid_undulation and Height_anomaly alone can be computed simultaneously.');
                end
            end
            if GUI==0 && (~isnumeric(volbapar) || any(rem(volbapar,1)~=0) || any(volbapar<2) || any(volbapar>25))
                error('The variable "Functional" must be an integer from the interval: 2, ..., 25.')
            end

            if volbagridcheck==1
                if coord==1 %Entered spherical coordinates                                                
                    if any(volbapar==10) || any(volbapar==23) 
                        if GUI==1 %If working with the GUI
                            errordlg('Geoid_undulation and Height_anomaly cannot be computed in the spherical coordinates. Reference surface to which geoid undulation and height anomaly are expressed is conventionally the surface of the reference ellipsoid. Grid that has been entered in the spherical coordinates (constant value of the radius r) does not refer to this surface. If you wish to compute these functionals, please use ellipsoidal type of the input coordinates and set value in the array Height above the reference surface (m) to zero.','Calculated parameters and output selection')
                        end
                        error('Geoid_undulation and Height_anomaly cannot be computed in the spherical coordinates. Reference surface to which geoid undulation and height anomaly are expressed is conventionally the surface of the reference ellipsoid. Grid that has been entered in the spherical coordinates (constant value of the radius r) does not refer to this surface. If you wish to compute these functionals, please use ellipsoidal type of the input coordinates and set value in the array Height above the reference surface (m) to zero.')
                    end
                end
            elseif volbadiskcheck==1
                if coord==1 %Entered spherical coordinates                                                
                    if any(volbapar==10) || any(volbapar==23)
                        if GUI==1 %If working with the GUI
                            errordlg('Geoid_undulation and Height_anomaly cannot be computed in the spherical coordinates. Reference surface to which geoid undulation and height anomaly are expressed is conventionally the surface of the reference ellipsoid. If you wish to compute these functionals, please use ellipsoidal type of the input coordinates and set value(s) in the array Ellipsoidal height/Spherical radius (m) to zero.','Calculated parameters and output selection')
                        end
                        error('Geoid_undulation and Height_anomaly cannot be computed in the spherical coordinates. Reference surface to which geoid undulation and height anomaly are expressed is conventionally the surface of the reference ellipsoid. If you wish to compute these functionals, please use ellipsoidal type of the input coordinates and set value(s) in the array Ellipsoidal height/Spherical radius (m) to zero.')
                    end
                end
            end
            
            if volbapar1==1 && volbapar2==1 && volbapar3==1 && volbapar4==1
                if GUI==1 %If working with the GUI
                    errordlg('Please choose at least one functional of the geopotential.',...
                        'Error in calculated paramters and output selection')
                end
                error('Please choose at least one functional of the geopotential.')   
            end 

            %If it not specified, GrafLab will not display data
            if isempty(display_data)
                display_data=0;
            end
            
            %Specifying the output folder and file
            if isempty(outname)
                if GUI==1 %Working with the GUI
                    set(findobj('tag','hlasky'),'string','Select output folder and output file',...
                           'fontsize',8,'foregroundcolor','k'); drawnow;

                    warn2=warndlg('Output folder and output file were not specified. Click OK and then select an output folder and output file. After the selection, the computation will start.');
                    waitfor(warn2);

                    [outname,outadresar]=uiputfile('*.*');
                    if outname==0
                        if GUI==1 %If working with the GUI
                            errordlg('Output folder and output file must be specified!');
                        end

                        set(findobj('tag','hlasky'),'string','',...
                           'fontsize',8,'foregroundcolor','k'); drawnow;

                        error('Output folder and output file must be specified!');
                    else
                        if find(outname=='.')>0
                            outname=outname(1:(find(outname=='.', 1,'last')-1));
                        end
                    end
                    
                    set(findobj('tag','OutFile'),'string',outname); %Display name of the output file
                    
                    set(findobj('tag','hlasky'),'string','',...
                           'fontsize',8,'foregroundcolor','k'); drawnow;
                else %Working without the GUI
                    error('The variable "Output_path" is empty.');
                end
            end
            
            %Some error checks of the "Export_data_txt" variable
            if isempty(Export_data_txt) || length(Export_data_txt)>1 || ~isnumeric(Export_data_txt) || rem(Export_data_txt,1)~=0
                error('The variable "Export_data_txt" must be an integer (0 -- do not export data into a txt file, 1 -- export data into a txt file).')
            end
            
            if Export_data_txt==1 || Export_data_txt==0
            else
                error('The variable "Export_data_txt" must be an integer (0 -- do not export data into a txt file, 1 -- export data into a txt file).')
            end
            
            %Some error checks of the "Export_report" variable
            if isempty(Export_report) || length(Export_report)>1 || ~isnumeric(Export_report) || rem(Export_report,1)~=0
                error('The variable "Export_report" must be an integer (0 -- do not export report, 1 -- export report).')
            end
            
            if Export_report==1 || Export_report==0
            else
                error('The variable "Export_report" must be an integer (0 -- do not export report, 1 -- export report).')
            end
            
            %Some error checks of the "Export_data_mat" variable
            if isempty(Export_data_mat) || length(Export_data_mat)>1 || ~isnumeric(Export_data_mat) || rem(Export_data_mat,1)~=0
                error('The variable "Export_data_mat" must be an integer (0 -- do not export data into a mat file, 1 -- export data into a mat file).')
            end
            
            if Export_data_mat==1 || Export_data_mat==0
            else
                error('The variable "Export_data_mat" must be an integer (0 -- do not export data into a mat file, 1 -- export data into a mat file).')
            end

            %Check of the entered values of DPI and number of colors
            if display_data==1 || display_data==2
                if isnan(DPI)==1 || DPI<0
                    if GUI==1 %If working with the GUI
                        errordlg('The value of "DPI" must be larger than zero.',...
                            'Display data settings');
                    end
                    error('The value of "DPI" must be larger than zero.');
                elseif isnan(ncolor)==1 || ncolor<2
                    if GUI==1 %If working with the GUI
                        errordlg('The number of colors must be larger than 1.',....
                            'Display data settings');
                    end
                    error('The number of colors must be larger than 1.');
                end
            end
            
            %Check if the computation of tensors in the LNOF using the MFCM 
            %has been selected (not allowed)
            if any(volbapar==8) || any(volbapar==9) || any(volbapar==14) || any(volbapar==15)
                if volbaALFs==2 && volbagridcheck==1 && STD==0
                    if GUI==1 %If working with the GUI
                        errordlg('The following functionals of the geopotential cannot be computed on a grid using the modified forward column method combined with Horner''s scheme: Disturbing_tensor_Txx_Tyy_Tzz, Disturbing_tensor_Txy_Txz_Tyz, Gravitational_tensor_Vxx_Vyy_Vzz and Gravitational_tensor_Vxy_Vxz_Vyz. In the case of a high maximum degree (1800 and higher), it is recommended to use the extended-range arithmetic approach. However, for the point-wise computation of the above mentioned functionals, the modified forward column method combined with Horner''s scheme can be applied.',...
                            'Calculated parameters and output selection');
                    end
                    error('The following functionals of the geopotential cannot be computed on a grid using the modified forward column method combined with Horner''s scheme: Disturbing_tensor_Txx_Tyy_Tzz, Disturbing_tensor_Txy_Txz_Tyz, Gravitational_tensor_Vxx_Vyy_Vzz and Gravitational_tensor_Vxy_Vxz_Vyz. In the case of a high maximum degree (1800 and higher), it is recommended to use the extended-range arithmetic approach. However, for the point-wise computation of the above mentioned functionals, the modified forward column method combined with Horner''s scheme can be applied.');
                end
            end

            %Check if display data in point-wise or load data approach has been
            %selected (not allowed)
            if (display_data==1 || display_data==2) && volbagridcheck~=1
                if GUI==1 %If working with the GUI
                    warn4=warndlg('Only data computed on a grid can be displayed. After clicking OK, the computation will start, but the computed data will not be displayed.');
                    waitfor(warn4);
                else
                    warning('Only data computed on a grid can be displayed. The computation will continue, but the computed data will not be displayed.'); %#ok<WNTAG>
                end
                display_data=0;
            end           
            
            if STD==0
              
                tic %Start clock to measure computation time

                %Loading GGM
                if isempty(GGMname) %Error message, if GGM file has not been imported
                    if GUI==1 %If working with the GUI
                        errordlg('Please input geopotential model file.',...
                            'Error in point type selection');
                    end
                    error('Please input geopotential model file.')
                end

                if GUI==1 %If working with the GUI
                    set(findobj('tag','hlasky'),'string',...
                        'Loading GGM file...','fontsize',8,...
                        'foregroundcolor','k'); drawnow;
                end
                
                if GUI==0 && Status_bar==1 %If working without the GUI
                    fprintf('Loading GGM file...\n')
                end
                
                if GUI==0
                    if strcmp(nmax,'nmaxGGM')
                        use_nmax_GGM=1;
                    else
                        use_nmax_GGM=0;
                    end
                end

                %Loading GGM
                if ~exist([GGMadresar,GGMname],'file') %Check whether the input GGM file exists
                    if GUI==1 %Working with the GUI
                        errordlg('The entered global geopotential model file does not exist.',...
                            'Error in global geopotential model and reference system selection');
                    end
                    error('The entered global geopotential model file does not exist.')
                end
                if strcmp(GGMname(end-3:end),'.gfc') %Input data in ICGEM format 
                    fGGMid=fopen([GGMadresar,GGMname]);

                    GM_old=GM;
                    R_old=R;
                    
                    cont=true;
                    while(~feof(fGGMid))                       
                        s=fgetl(fGGMid); 
                        if ~strncmpi(s,'product_type',12) && cont==true
                            continue
                        else
                            cont=false;
                        end
                        if strncmpi(s,'earth_gravity_constant',22)
                            s_temp=isspace(s);
                            s=s(find(s_temp>0, 1 ):end);
                            s=strtrim(s);
                            s_temp=isspace(s);
                            if any(s_temp>0)
                                s=s(1:find(s_temp==1, 1 ));
                            end
                            GM=str2num(s); %#ok<*ST2NM>
                        end
                        if strncmpi(s,'radius',6)
                            s_temp=isspace(s);
                            s=s(find(s_temp>0, 1 ):end);
                            s=strtrim(s);
                            s_temp=isspace(s);
                            if any(s_temp>0)
                                s=s(1:find(s_temp==1, 1 ));
                            end
                            R=str2num(s); %#ok<*ST2NM>
                        end
                        if strncmpi(s,'norm',4)
                            s_temp=isspace(s);
                            s=s(find(s_temp>0, 1 ):end);
                            s=strtrim(s);
                            if ~strncmpi(s,'fully_normalized',16)
                                if GUI==1 %If working with the GUI
                                    errordlg('GrafLab can work with *.gfc files only if the coefficients are fully normalized.',...
                                        'Geopotential model and reference system selection');
                                end
                                error('GrafLab can work with *.gfc files only if the coefficients are fully normalized.')
                            end
                        end
                        if strncmpi(s,'end_of_head',11)
                            break
                        end
                    end
                    clear s s_temp
                    
                    if GM~=GM_old || R~=R_old
                        if GUI==1 %If working with the GUI
                            warn5=warndlg(sprintf('The value of the geocentric gravitational constant (GM) and/or the radius of the reference sphere (R) in the header of the input "%s" file are different from those entered through the GUI. After clicking the "OK" button, GrafLab will automatically use the GM and R values from the "%s" file and will continue in the computation.',GGMname,GGMname));
                            waitfor(warn5);
                        elseif GUI==0 %If working without the GUI
                            warning('The value of the geocentric gravitational constant (GM) and/or the radius of the reference sphere (R) in the header of the input "%s" file are different from those entered through the "GrafLab.m" function. GrafLab will now automatically use the GM and R values from the "%s" file and will continue in the computation.',GGMname,GGMname); %#ok<WNTAG>
                        end
                    end
                    
                    GGM=textscan(fGGMid,'%s%f%f%f%f%f%f');

                    if any(strcmp(GGM{1},'gfc')==0)
                        if GUI==1 %If working with the GUI
                            errordlg('Unsupported format of the GGM file.',...
                                'Geopotential model and reference system selection');
                        end
                        error('Unsupported format of the GGM file.')
                    end
                    
                    fclose(fGGMid);

                    GGM=GGM(2:end);
                    GGM=cell2mat(GGM);
                    GGM_length=1;
                elseif strcmp(GGMname(end-3:end),'.mat') %Input data in MAT format 
                    GGM=load([GGMadresar,GGMname]);
                    
                    GGM_length=numel(fieldnames(GGM));
                    if GGM_length==1
                        GGM=struct2cell(GGM);
                    elseif GGM_length==3
                        nminGGM=GGM.nmin;
                        GGM=rmfield(GGM,'nmin');
                        nmaxGGM=GGM.nmax;
                        GGM=rmfield(GGM,'nmax');
                        GGM=struct2cell(GGM);
                        
                    else
                        if GUI==1 %If working with the GUI
                        errordlg('Unsupported format of the GGM file.',...
                            'Geopotential model and reference system selection');
                        end
                        error('Unsupported format of the GGM file.')
                    end
                    GGM=cell2mat(GGM);
                else
                    GGM=load([GGMadresar,GGMname]);
                    GGM_length=1;
                end
                
                if GGM_length==1
                    [rows_GGM,cols_GGM]=size(GGM);
                    if cols_GGM<4
                        if GUI==1 %If working with the GUI
                            errordlg('Unsupported format of the GGM file.',...
                                'Geopotential model and reference system selection');
                        end
                        error('Unsupported format of the GGM file.')
                    end
                    
                    GGM=GGM(:,1:4);
                    GGM=sortrows(GGM,1);
                    stupen=GGM(:,1);
                    rad=GGM(:,2);
                    C=GGM(:,3);
                    S=GGM(:,4);
                    nminGGM=min(stupen);
                       
                    if max(stupen)<2
                        if GUI==1 %If working with the GUI
                            errordlg('The maximum degree of the improted GGM file must be at least 2.',...
                                'Geopotential model and reference system selection');
                        end
                        error('The maximum degree of the improted GGM file must be at least 2.')
                    end
                    
                    if stupen(1)==0 && stupen(2)==2
                        %Some geopotential models include the zero-degree
                        %term, but omit coefficients of degree 1 if they
                        %are zero. If this is the case, GrafLab add the
                        %zero coefficients of degree one into the model
                        del10=find(stupen==1 & rad==0,1); %If the coefficients of the degree 1 and order 0 are missing
                        if isempty(del10)
                            stupen=[stupen(1);1;stupen(2:end)];
                            rad=[rad(1);0;rad(2:end)];
                            C=[C(1);0;C(2:end)];
                            S=[S(1);0;S(2:end)];
                        end                             

                        del11=find(stupen==1 & rad==1,1); %If the coefficients of the degree 1 and order 1 are missing             
                        if isempty(del11)
                            stupen=[stupen(1:2);1;stupen(3:end)];
                            rad=[rad(1:2);1;rad(3:end)];
                            C=[C(1:2);0;C(3:end)];
                            S=[S(1:2);0;S(3:end)];
                        end
                        clear del00 del10 del11
                    end
                    
                    nmaxGGM=max(stupen);
                    
                    %Identification of GGM file format
                    if nminGGM==0
                        if stupen(1)==0 && stupen(2)==1 && stupen(3)==1 && stupen(4)==2 && stupen(5)==2 && rad(1)==0 && rad(2)==0 && rad(3)==1 && rad(4)==0 && rad(5)==1
                        else
                            if GUI==1 %If working with the GUI
                                errordlg('Wrong format of the input GGM file.',...
                                    'Geopotential model and reference system selection');
                            end
                            error('Wrong format of the input GGM file.')
                        end
                    end
                elseif GGM_length==3
                    [rows_GGM,cols_GGM]=size(GGM);
                    if cols_GGM<2
                        if GUI==1 %If working with the GUI
                            errordlg('Unsupported format of the GGM file.',...
                                'Geopotential model and reference system selection');
                        end
                        error('Unsupported format of the GGM file.')
                    end
                    
                    C=GGM(:,1);
                    S=GGM(:,2);
                    
                    sum_nmin=0;
                    for n=0:nminGGM
                        sum_nmin=sum_nmin+n;
                    end
                    sum_nmax=0;
                    for n=0:(nmaxGGM+1)
                        sum_nmax=sum_nmax+n;
                    end
                    stupen=zeros(sum_nmax-sum_nmin,1);
                    rad=stupen;
                    nn=1;
                    for n=nminGGM:nmaxGGM
                        for m=0:n
                            stupen(nn)=n;
                            rad(nn)=m;
                            nn=nn+1;
                        end
                    end
                end
                clear GGM
                
                set(findobj('tag','hlasky'),'string',...
                        '','fontsize',8,'foregroundcolor','k'); drawnow;
                
                %Error messages for the nmin value
                if isempty(nmin) || length(nmin)>1 || ~isnumeric(nmin) || rem(nmin,1)~=0
                    if GUI==1 %Working with the GUI
                        errordlg('The value of "nmin" must be an integer.',...
                            'Error in geopotential model and reference system selection')
                    end
                    error('The value of "nmin" must be an integer.')
                end
                if nmin<0
                    if GUI==1 %Working with the GUI
                        errordlg('The value of "nmin" cannot be negative.',...
                            'Error in geopotential model and reference system selection')
                    end
                    error('The value of "nmin" cannot be negative.')
                elseif nmin>nmaxGGM
                    if GUI==1 %Working with the GUI
                        errordlg('The value of "nmin" exceedes the "nmax" value of GGM.',...
                            'Error in geopotential model and reference system selection')
                    end
                    error('The value of "nmin" exceedes the "nmax" value of GGM.')
                end
                if isnan(nmin)==1
                    if GUI==1 %Working with the GUI
                        errordlg('Please input the "nmin" value.',...
                            'Error in geopotential model and reference system selection')
                    end
                    error('Please input the "nmin" value.')
                end
                
                %Value of nmax and error messages   
                if use_nmax_GGM==1 %If GrafLab takes automatically the maximum degree for the 
                    %synthesis from the global geopotential model file
                    nmax=nmaxGGM;
                end
                if isempty(nmax) || length(nmax)>1 || ~isnumeric(nmax) || rem(nmax,1)~=0
                    if GUI==1 %Working with the GUI
                        errordlg('The value of "nmax" must be an integer.',...
                            'Error in geopotential model and reference system selection')
                    end
                    error('The value of "nmax" must be an integer.')
                end
                if nmax>nmaxGGM
                    if GUI==1 %Working with the GUI
                        errordlg('The entered value of "nmax" exceedes the "nmax" value of GGM.',...
                            'Error in geopotential model and reference system selection')
                    end
                    error('The entered value of "nmax" exceedes the "nmax" value of GGM.')
                elseif nmax<2
                    if GUI==1 %Working with the GUI
                        errordlg('The value of "nmax" must be at least 2.',...
                            'Error in geopotential model and reference system selection')
                    end
                    error('The value of "nmax" must be at least 2.')
                elseif nmin>nmax
                    if GUI==1 %Working with the GUI
                        errordlg('The value of "nmin" cannot be larger than "nmax" value.',...
                            'Error in geopotential model and reference system selection')
                    end
                    error('The value of "nmin" cannot be larger than "nmax" value.')
                end
                if isnan(nmax)==1
                    if GUI==1 %Working with the GUI
                        errordlg('Please input the "nmax" value.',...
                            'Error in geopotential model and reference system selection')
                    end
                    error('Please input the "nmax" value.')
                end
                
                if nmin>0
                    if any(volbapar==20)
                        if GUI==1 %If working with the GUI
                            errordlg('Gravity disturbance cannot be computed if "nmin>0".',...
                                'Calculated parameters and output selection');
                        end
                        error('Gravity disturbance cannot be computed if "nmin>0".');
                    end
                end
                
                if nmin>0
                    if any(volbapar==10) || any(volbapar==23)
                        if GUI==1 %If working with the GUI
                            errordlg('Geoid_undulation and Height_anomaly cannot be computed if "nmin>0".');
                        end
                        error('Geoid_undulation and Height_anomaly cannot be computed if "nmin>0".');
                    end
                    
                    if any(volbapar==9) || any(volbapar==15)
                        if GUI==1 %If working with the GUI
                            errordlg('The following functionals of the geopotential cannot be computed if "nmin>0": Disturbing_tensor_Txy_Txz_Tyz and Gravitational_tensor_Vxy_Vxz_Vyz.',...
                                'Calculated parameters and output selection');
                        end
                        error('The following functionals of the geopotential cannot be computed if "nmin>0": Disturbing_tensor_Txy_Txz_Tyz and Gravitational_tensor_Vxy_Vxz_Vyz.');
                    end
                    
                    if any(volbapar==8) || any(volbapar==14) 
                        LNOFnmin=1; %Logical 1
                        
                        if any(volbapar~=8 & volbapar~=9 & volbapar~=14 & volbapar~=15 & volbapar~=1)                                              
                            if GUI==1 %If working with the GUI
                                errordlg('The following functionals of the geopotential cannot be computed simultaneously with other functionals if "nmin>0": Disturbing_tensor_Txx_Tyy_Tzz and Gravitational_tensor_Vxx_Vyy_Vzz.',...
                                    'Calculated parameters and output selection');
                            end
                            error('The following functionals of the geopotential cannot be computed simultaneously with other functionals if "nmin>0": Disturbing_tensor_Txx_Tyy_Tzz and Gravitational_tensor_Vxx_Vyy_Vzz.');
                        end
                    else
                        LNOFnmin=0; %Logical 0
                    end
                end
                
                if nminGGM>0 && nminGGM~=nmin && (any(volbapar==8) || any(volbapar==14))
                    if GUI==1 %If working with the GUI
                        errordlg('The functionals Disturbing_tensor_Txx_Tyy_Tzz and Gravitational_tensor_Vxx_Vyy_Vzz cannot be computed if the minimum degree of the imported GGM file is larger than zero and at the same time the minimum degree of the GGM is different from the one that is to be used in the synthesis (specified in GrafLab).',...
                            'Calculated parameters and output selection');
                    end
                    error('The functionals Disturbing_tensor_Txx_Tyy_Tzz and Gravitational_tensor_Vxx_Vyy_Vzz cannot be computed if the minimum degree of the imported GGM file is larger than zero and at the same time the minimum degree of the GGM is different from the one that is to be used in the synthesis (specified in GrafLab).');
                end
                
                if nmax<nmaxGGM %If the maximum degree that is to be used in the
                    %synthesis is smaller than the maximum degree of the
                    %imported GGM file, the coefficients beyond the former
                    %can be deleted to save some RAM
                    idx_del=stupen>nmax;
                    stupen(idx_del)=[];
                    rad(idx_del)=[];
                    C(idx_del)=[];
                    S(idx_del)=[];
                    clear idx_del
                end
                
                %Loading DMR
                if any(volbapar==10) || any(volbapar==23)
                                       
                    if GUI==1 %If working with the GUI
                        set(findobj('tag','hlasky'),'string','Please select DTM file',...
                            'fontsize',8,'foregroundcolor','k'); drawnow;
                        [loadnameDMR,loadadresarDMR]=uigetfile('*.*','Select File Containg Spherical Harmonic Coefficients of the Topography');
                    end

                    if loadnameDMR==0
                        if GUI==1 %If working with the GUI
                            errordlg('To compute Geoid_undulation/Height_anomaly, DTM file must be imported!',...
                                'Error in geopotential model and reference system selection');
                        end
                        error('To compute Geoid_undulation/Height_anomaly, DTM file must be imported!')
                    else
                        
                        if GUI==1 %If working with the GUI
                            set(findobj('tag','hlasky'),'string',...
                                'Loading DTM file...','fontsize',8,...
                                'foregroundcolor','k'); drawnow;
                        end
                        
                        if GUI==0 && Status_bar==1 %If working without the GUI
                            fprintf('Loading DTM file...\n');
                            
                            if isempty(DTM_path)
                                error('The variable "DTM_path" cannot be empty if computing functional "Geoid_undulation" or "Height_anomaly."')
                            end
                        end

                        if ~exist([loadadresarDMR,loadnameDMR],'file') %Check whether the input DTM file exists
                            if GUI==1 %Working with the GUI
                                errordlg('The entered DTM file does not exist.',...
                                    'Error in global geopotential model and reference system selection');
                            end
                            error('The entered DTM file does not exist.')
                        end
                        if strcmp(loadnameDMR(end-3:end),'.mat')
                            DMR=load([loadadresarDMR,loadnameDMR]);
                            DMR=struct2cell(DMR);
                            DMR=cell2mat(DMR);                           
                        else
                            DMR=load([loadadresarDMR,loadnameDMR]);
                        end
                       
                        [rows_DMR,cols_DMR]=size(DMR); %#ok<*ASGLU>
                        if cols_DMR<4
                            if GUI==1 %If working with the GUI
                                errordlg('Wrong format of the input DTM file.',...
                                    'Geopotential model and reference system selection');
                            end
                            error('Wrong format of the input DTM file.')
                        end
                       
                        DMR=DMR(:,1:4);
                        
                        if max(DMR(:,1))<nmax
                            if GUI==1 %If working with the GUI
                                errordlg('The maximum degree in the DTM file is smaller than the used "nmax" value of the geopotential model.',...
                                    'Geopotential model and reference system selection');
                            end
                            error('The maximum degree in the DTM file is smaller than the used "nmax" value of the geopotential model.')
                        end
                        
                        if GUI==1 %If working with the GUI
                            set(findobj('tag','hlasky'),'string','',...
                                'foregroundcolor','k'); drawnow;
                        end
                    end
                end  
                
                %Sorting spherical harmonic coeffcients of DTM
                if any(volbapar==10) || any(volbapar==23)
                    DMR(DMR(:,1)>nmaxGGM,:)=[];

                    if DMR(1,1)==0 && DMR(2,1)==1 && DMR(3,1)==1 && DMR(4,1)==2 && DMR(1,2)==0 && DMR(2,2)==0 && DMR(3,2)==1 && DMR(4,2)==0
                    else
                        DMR=sortrows(DMR,1);
                    end
                        
                    HC=DMR(:,3);
                    HS=DMR(:,4);

                    clear DMR
                end
                
                %Indices of the spherical harmonic coefficients
                index=zeros(nmax+1-nminGGM,1);
                index(1)=1;
                for i=1:(nmax-nminGGM)
                    index(i+1)=index(i)+nminGGM+i;
                end
                
                
                %% Computation of functionals of the geopotential on a regular grid
                if volbagridcheck==1

                    %Entered coordinates of the grid
                    if GUI==1 %If working with the GUI
                        fimin=str2num(get(findobj('tag','fimin'),'string'));
                        fistep=str2num(get(findobj('tag','fistep'),'string'));
                        fimax=str2num(get(findobj('tag','fimax'),'string'));
                        lambdamin=str2num(get(findobj('tag','lambdamin'),'string'));
                        lambdastep=str2num(get(findobj('tag','lambdastep'),'string'));
                        lambdamax=str2num(get(findobj('tag','lambdamax'),'string'));
                        h=str2num(get(findobj('tag','hgrid'),'string'));
                    elseif GUI==0 %If working without the GUI
                        fimin=lat_min; 
                        fistep=lat_step;
                        fimax=lat_max;
                        lambdamin=lon_min; 
                        lambdastep=lon_step;
                        lambdamax=lon_max;
                        clear lat_min lat_step lat_max lon_min lon_step lon_max
                    end

                    %Input coordinates and some error checks
                    if GUI==0 && strcmp(fistep,'empty') && strcmp(fimax,'empty')
                        
                        if isempty(fimin) || ~isnumeric(fimin) || ~isvector(fimin)
                            error('When using a vector to define latitudes, the variable "lat_min" must be a real-valued scalar of vector.')
                        end
                        
                        %Vector of latitudes defined via the variable
                        %"lat_min"
                        fi=fimin(:); 
                        
                    else
                        if isempty(fimin) || ~isnumeric(fimin) || length(fimin)>1
                            if GUI==1
                                errordlg('The value of "Lat. min" must be a real-valued scalar.',...
                                    'Error in irregular surface selection panel');
                                error('The value of "Lat. min" must be a real-valued scalar.');
                            elseif GUI==0
                                error('The variable "lat_min" must be a real-valued scalar.');
                            end
                        end
                        if isempty(fistep) || ~isnumeric(fistep) || length(fistep)>1
                            if GUI==1
                                errordlg('The value of "Lat. step" must be a real-valued scalar.',...
                                    'Error in irregular surface selection panel');
                                error('The value of "Lat. step" must be a real-valued scalar.');
                            elseif GUI==0
                                error('The variable "lat_step" must be a real-valued scalar.');
                            end
                        end
                        if isempty(fimax) || ~isnumeric(fimax) || length(fimax)>1
                            if GUI==1
                                errordlg('The value of "Lat. max" must be a real-valued scalar.',...
                                    'Error in irregular surface selection panel');
                                error('The value of "Lat. max" must be a real-valued scalar.');
                            elseif GUI==0
                                error('The variable "lat_max" must be a real-valued scalar.');
                            end
                        end

                        if fimin>fimax
                            if GUI==1
                                errordlg('The value of "Lat. min" must be smaller than the "Lat. max" value.',...
                                    'Error in irregular surface selection panel');
                                error('The value of "Lat. min" must be smaller than the "Lat. max" value.'); 
                            elseif GUI==0
                                error('The value of "lat_min" must be smaller than the "lat_max" value.'); 
                            end
                        end
                        if fistep<=0
                            if GUI==1
                                errordlg('The value of "Lat. step" must be larger than zero.',...
                                    'Error in irregular surface selection panel');
                                error('The value of "Lat. step" must be larger than zero.');
                            elseif GUI==0
                                error('The value of "lat_step" must be larger than zero.');
                            end
                        end

                        if fimin>90 || fimin<-90
                            if GUI==1
                                errordlg('The value of "Lat. min" must be within the interval <-90 deg, 90 deg>.',...
                                    'Error in irregular surface selection panel');
                                error('The value of "Lat. min" must be within the interval <-90 deg, 90 deg>.');
                            elseif GUI==0
                                error('The value of "lat_min" must be within the interval <-90 deg, 90 deg>.');
                            end
                        end
                        if fimax>90 || fimax<-90
                            if GUI==1
                                errordlg('The value of "Lat. max" must be within the interval <-90 deg, 90 deg>.',...
                                    'Error in irregular surface selection panel');
                                error('The value of "Lat. max" must be within the interval <-90 deg, 90 deg>.');
                            elseif GUI==0
                                error('The value of "lat_max" must be within the interval <-90 deg, 90 deg>.');
                            end
                        end

                        fi=(fimin:fistep:fimax)';
                    end              
                    
                    if GUI==0 && strcmp(lambdastep,'empty') && strcmp(lambdamax,'empty')
                        
                        if isempty(lambdamin) || ~isnumeric(lambdamin) || ~isvector(lambdamin)
                            error('When using a vector to define longitudes, the variable "lon_min" must be a real-valued scalar of vector.')
                        end
                        
                        %Vector of latitudes defined via the variable
                        %"lat_min"
                        lambda=lambdamin(:); 
                        
                    else
                        if isempty(lambdamin) || ~isnumeric(lambdamin) || length(lambdamin)>1
                            if GUI==1
                                errordlg('The value of "Lon. min" must be a real-valued scalar.',...
                                    'Error in irregular surface selection panel');
                                error('The value of "Lon. min" must be a real-valued scalar.');
                            elseif GUI==0
                                error('The variable "lon_min" must be a real-valued scalar.');
                            end
                        end
                        if isempty(lambdastep) || ~isnumeric(lambdastep) || length(lambdastep)>1
                            if GUI==1
                                errordlg('The value of "Lon. step" must be a real-valued scalar.',...
                                    'Error in irregular surface selection panel');
                                error('The value of "Lon. step" must be a real-valued scalar.');
                            elseif GUI==0
                                error('The variable "lon_step" must be a real-valued scalar.');
                            end
                        end
                        if isempty(lambdamax) || ~isnumeric(lambdamax) || length(lambdamax)>1
                            if GUI==1
                                errordlg('The value of "Lon. max" must be a real-valued scalar.',...
                                    'Error in irregular surface selection panel');
                                error('The value of "Lon. max" must be a real-valued scalar.');
                            elseif GUI==0
                                error('The variable "lon_max" must be a real-valued scalar.');
                            end
                        end

                        if lambdamin>lambdamax
                            if GUI==1
                                errordlg('The value of "Lon. min" must be smaller than Lon. max value.',...
                                    'Error in irregular surface selection panel');
                                error('The value of "Lon. min" must be smaller than Lon. max value.');
                            elseif GUI==0
                                error('The value of "lon_min" must be smaller than Lon. max value.');
                            end
                        end
                        if lambdastep<=0
                            if GUI==1
                                errordlg('The value of "Lon. step" must be larger than zero.',...
                                    'Error in irregular surface selection panel');
                                error('The value of "Lon. step" must be larger than zero.');
                            elseif GUI==0
                                error('The value of "lon_step" must be larger than zero.');
                            end
                        end

                        if lambdamin>360 || lambdamin<-180
                            if GUI==1
                                errordlg('The value of "Lon. min" must be within the interval <-180 deg, 180 deg> or <0 deg, 360 deg>.',...
                                    'Error in irregular surface selection panel');
                                error('The value of "Lon. min" must be within the interval <-180 deg, 180 deg> or <0 deg, 360 deg>.');
                            elseif GUI==0
                                error('The value of "lon_min" must be within the interval <-180 deg, 180 deg> or <0 deg, 360 deg>.');
                            end
                        end
                        if lambdamax>360 || lambdamax<-180
                            if GUI==1
                                errordlg('The value of "Lon. max" must be within the interval <-180 deg, 180 deg> or <0 deg, 360 deg>.',...
                                    'Error in irregular surface selection panel');
                                error('The value of "Lon. max" must be within the interval <-180 deg, 180 deg> or <0 deg, 360 deg>.');
                            elseif GUI==0
                                error('The value of "lon_max" must be within the interval <-180 deg, 180 deg> or <0 deg, 360 deg>.');
                            end
                        end
                        if (lambdamax-lambdamin)>360
                            if GUI==1
                                errordlg('The longitudes must be in the range <-180 deg, 180 deg> or <0 deg, 360 deg>.',...
                                    'Error in irregular surface selection panel');
                                error('The longitudes must be in the range <-180 deg, 180 deg> or <0 deg, 360 deg>.');
                            elseif GUI==0
                                error('The longitudes must be in the range <-180 deg, 180 deg> or <0 deg, 360 deg>.');
                            end
                        end
                        
                        lambda=(lambdamin:lambdastep:lambdamax)'; 
                    end 
                    
                    if isempty(h) || ~isnumeric(h) || length(h)>1
                        if GUI==1
                            errordlg('The value of "Height above the reference surface" must be a real-valued scalar.',...
                                'Error in irregular surface selection panel');
                            error('The value of "Height above the reference surface" must be a real-valued scalar.');
                        elseif GUI==0
                            error('The variable "h" must be a real-valued scalar.');
                        end
                    end
                    
                    length_fi=length(fi);
                    length_lambda=length(lambda);
                    %Is the input grid symetric with respect to the
                    %equator? If yes, GrafLab exploits the symetry property
                    %of Legendre functions which results in an increased
                    %computational speed.
                    if rem(length_fi,2)==0 %The number of latitudes
                            %is an even number
                        cond1=true;
                    else %The number of latitudes is an odd number
                        mid_lat_idx=ceil(length_fi/2);
                        if abs(fi(mid_lat_idx))<100*eps %The middle latitude
                            %is zero up to a given level of accuracy (here 100*eps)
                            cond1=true;
                        else
                            cond1=false;
                        end
                    end
                    %Next, GrafLab checks the symetry of the
                    %latitudes with respect to the equator (up to the given
                    %level of accuracy 100*eps)
                    if rem(length_fi,2)==0 %The number of latitudes
                            %is an even number
                        odd=0;
                        if max(abs(abs(fi(1:(length_fi/2)))-flipud(fi((length_fi/2+1):end))))<100*eps
                            cond2=true;
                        else
                            cond2=false;
                        end
                    else %The number of latitudes is an odd number
                        odd=1;
                        if max(abs(abs(fi(1:(mid_lat_idx-1)))-flipud(fi((mid_lat_idx+1):end))))<100*eps  
                            cond2=true;
                        else
                            cond2=false;
                        end
                    end
                    if cond1 && cond2
                        symmetric_grid=true; %The grid is symetric with respect to the equator
                    else
                        symmetric_grid=false; %The grid is not symetric with respect to the equator
                    end
                    
                    fi=pi/180*(fi(:));
                    lambda=pi/180*(lambda(:));

                    %Grid, which is to be displayed has to have at least two
                    %points in latitude parallels and at least two points in
                    %longitude parallels.
                    if display_data==1 || display_data==2
                        if length_fi<2 || length_lambda<2
                            if GUI==1 %If working with the GUI
                                warn3=warndlg('To display computed data on a grid, it must contain at least two distinct points in one parallel and two distinct points in one meridian. After clicking OK, the computation will start, but the data will not be displayed.');
                                waitfor(warn3);
                            elseif GUI==0 %If working without the GUI
                                warning('To display computed data on a grid, it must contain at least two distinct points in one parallel and two distinct points in one meridian. The computation will continue, but the data will not be displayed.'); %#ok<WNTAG>
                            end
                            display_data=0;
                        end
                    end
                    
                    if coord==1 %Entered spherical coordinates                       
                        %Spherical radius
                        r=(R+h)*ones(length_fi,1);
                        hsph=h;
                        
                        %Spherical latitude
                        fiG=fi;  
                        
                        %Transformation of spherical latitude into the ellipsoidal
                        %latitude
                        [X,Y,Z]=sph2cart(0*fiG,fiG,r);
                        [fi,lambda_del,h]=cart2ell(X,Y,Z,[aEl eEl]);

                        clear X Y Z lambda_del
                    elseif coord==0 %Entered ellipsoidal coordinates
                        %Trasformation of (fi, lambda, h) into (X, Y, Z)
                        [X,Y,Z]=ell2cart(fi,0*zeros(length_fi,1),h*ones(length_fi,1),[aEl eEl]);  
                        r=sqrt(X.*X+Y.*Y+Z.*Z); %Radius

                        %Spherical latitude
                        fiG=atan(Z./sqrt(X.*X+Y.*Y)); 
                        
                        clear X Y Z 
                    end

                    %Computation of the coefficients C0,0; C2,0; ...; C20,0
                    %of the selected ellipsoid
                    CEl=zeros(length(C),1);
                    for n=0:10
                        CEl(2*n==stupen & rad==0,1)=((-1)^n*(3*eEl^(2*n))/((2*n+1)*(2*n+3)*sqrt(4*n+1))*(1-n-5^(3/2)*n*CEl_20/eEl^2)).*(aEl./R).^(2*n).*(GMEl/GM);
                    end                                

                    if any(volbapar==11) || any(volbapar==12) || any(volbapar==13) || any(volbapar==14) || any(volbapar==15) || any(volbapar==16) || any(volbapar==17) || any(volbapar==18) || any(volbapar==20) || any(volbapar==25)
                        grav=1;
                    else
                        grav=0;
                    end

                    if  any(volbapar==2) || any(volbapar==3) || any(volbapar==4) || any(volbapar==5) || any(volbapar==6) || any(volbapar==7) || any(volbapar==8) || any(volbapar==9) || any(volbapar==10) || any(volbapar==19) || any(volbapar==21) || any(volbapar==22) || any(volbapar==23) || any(volbapar==24)
                        por=1;
                        deltaC=C-CEl;
                    else
                        por=0;
                    end
 
                    if any(volbapar==20)
                        normal=1;
                    else
                        normal=0;
                    end

                    clear stupen rad                   
                    if normal==0
                        clear CEl
                    end
                    
                    %Initialization
                    eta=0; ksi=0; Theta=0; T=0; T_rr=0; Trr=0; Trf=0; Trl=0; Tff=0; %#ok<*NASGU>
                    Tfl=0; Tll=0; Tzz=0; Txx=0; Tyy=0; N=0; V=0; Vrr=0; Vrf=0; Vrl=0; Vff=0;
                    Vfl=0; Vll=0; g=0; g_sa=0; W=0; anomalia_sa=0; porucha=0;
                    porucha_sa=0; zetaEl=0; zeta=0; Wrr=0; Wr=0; Wfi=0; Wlambda=0;
                    Ur=0; Ufi=0; N1c=0; N2c=0; H=0;

                    %If geoid/quasigeoid is to be computed
                    if any(volbapar==10) || any(volbapar==23)
                        geoid=1;
                        if h~=0
                            if GUI==1 %If working with the GUI
                                errordlg('To compute Geoid_undulation or Height_anomaly, the ellipsoidal height must be set to zero.',...
                                    'Error in point type selection');
                            end
                            error('To compute Geoid_undulation or Height_anomaly, the ellipsoidal height must be set to zero.');
                        end
                    else
                        geoid=0;
                    end

                    %Initialization of the matrices and vectors for the computation of fnALFs
                    length_fiG=length(fiG);
                    if symmetric_grid==true && volbaALFs==3 %GrafLab takes advantage
                        %of the symmetric grid only when extended range
                        %arithmetic is used to compute Legednre functions
                        length_fiG_ALFs=ceil(length_fiG/2);
                        length_fi_ALFs=ceil(length_fi/2);
                    else
                        length_fiG_ALFs=length_fiG;
                        length_fi_ALFs=length_fi;
                    end
                    Pnm=zeros(length_fiG,nmax+1);
                    q_all=R./r;
                    q=q_all(1:length_fiG_ALFs);
                    q2=q.^2;
                    u_all=cos(fiG);
                    u=u_all(1:length_fiG_ALFs);
                    t_all=sin(fiG);
                    t=t_all(1:length_fiG_ALFs);
                    
                    %Initialization for extended-range arithmetic approach
                    if volbaALFs==3
                                               
                        bit=mexext; %Bit version of Matlab
                        bit=bit(end-1:end);
                        bit=str2double(bit);
                        if bit==32
                            bit=32;
                        elseif bit==64
                            bit=64;
                        else
                            bit=64;
                        end
                        
                        nmax23=nmax*2+3;
                        rr=zeros(nmax23,1); ri=rr;
                        dd=zeros(nmax,1); am=dd; bm=am;

                        m1=1:nmax23;
                        rr(m1)=sqrt(m1);
                        ri(m1)=1./rr;
                        m2=1:nmax;
                        dd(m2)=rr(2*m2+3).*ri(2*m2+2);

                        IND=960;
                        BIG=2^IND;
                        BIGI=2^(-IND);
                        BIGS=2^(IND/2);
                        BIGSI=2^(-IND/2);
                        ROOT3=1.732050807568877;
                        
                        if bit==32
                            pm=am;
                            ps1=zeros(length_fiG_ALFs,nmax); 
                            ips1=ps1;
                            x=ROOT3*u.*q;
                            ix=zeros(size(x));
                            ps1(:,1)=x;
                            ips1(:,1)=ix;
                            for m3=2:nmax
                                x=(dd(m3-1)*u).*x.*q;
                                y=abs(x);
                                iy=y>=BIGS;
                                if any(iy)
                                    x(iy)=x(iy)*BIGI;
                                    ix(iy)=ix(iy)+1;
                                end
                                iy=y<BIGSI;
                                if any(iy)
                                    x(iy)=x(iy)*BIG;
                                    ix(iy)=ix(iy)-1;
                                end
                                ps1(:,m3)=x;
                                ips1(:,m3)=ix;
                            end
                        elseif bit==64
                            tq=t.*q;
                            temp1=zeros(length_fiG_ALFs,1);
                            temp2=ones(length_fiG_ALFs,1);
                            temp3=temp2;
                            temp4=temp1;
                            temp5=temp1+BIGI;
                            ps1b=zeros(length_fiG_ALFs,nmax); 
                            ips1b=ps1b;
                            xb=ROOT3*u.*q;
                            ixb=zeros(size(xb));
                            ps1b(:,1)=xb;
                            ips1b(:,1)=ixb;
                            for m3=2:nmax
                                xb=(dd(m3-1)*u).*xb.*q;
                                yb=abs(xb);
                                iyb=yb>=BIGS;
                                if any(iyb)
                                    xb(iyb)=xb(iyb)*BIGI;
                                    ixb(iyb)=ixb(iyb)+1;
                                end
                                iyb=yb<BIGSI;
                                if any(iyb)
                                    xb(iyb)=xb(iyb)*BIG;
                                    ixb(iyb)=ixb(iyb)-1;
                                end
                                ps1b(:,m3)=xb;
                                ips1b(:,m3)=ixb;
                            end
                        end
                        
                        clear dd
                    end
                    u=u_all;
                    t=t_all;
                    
                    %Initialization of the matrices and vectors for the 
                    %computation of the first-order derivatives of fnALFs
                    if any(volbapar==3) || any(volbapar==4) || any(volbapar==6) || any(volbapar==7) || any(volbapar==12) || any(volbapar==13) || any(volbapar==16)  || any(volbapar==20)
                        dALFs=1;
                        dPnm=zeros(length_fiG,nmax+1);
                        qu=q_all./u;
                        tu=t./u;
                        
                        %Treatment of the dPnm singularity
                        singdPnm=fi==pi/2 | fi==-pi/2;
                    else
                        dALFs=0;
                    end   
                    
                    %Initialization of the matrices and vectors for the 
                    %computation of the second-order derivatives of fnALFs
                    if any(volbapar==6) || any(volbapar==12)
                        ddALFs=1;
                        ddPnm=zeros(length_fiG,nmax+1);
                        
                        %Treatment of the ddPnm singularity
                        singddPnm=fi==pi/2 | fi==-pi/2;
                    else
                        ddALFs=0;
                    end   
                    
                    %Status line
                    progressbar=findobj('tag','hlasky');

                    if GUI==0 && Status_bar==1 %If working without the GUI
                        fprintf('Progress: m = ')
                    end
                    
                    %% Summation over m
                    for m=nmax:-1:0

                        %Update of the progress bar
                        if GUI==1 && rem(m,10)==0 %If working with the GUI
                            set(progressbar,'string',...
                                sprintf('Progress: m = %5.0d',m),...
                                'fontsize',8); drawnow;
                        end
                        
                        if GUI==0 && Status_bar==1 %If working without the GUI
                            if rem(m,10)==0 || m==0
                                fprintf('%d,',m)
                            end
                        end
                        
                        m_min=1+max([0 nminGGM-m]);

                        %Selection of the spherical harmonic coefficients of order m
                        %======================================================
                        if grav==1 %C's spherical harmonic coefficients for the functionals without the normal gravity field
                            Cm=C(index((m+m_min-nminGGM):end)+m);
                        end

                        if por==1 %C's spherical harmonic coefficients for the functionals with the disturbing field
                            deltaCm=deltaC(index((m+m_min-nminGGM):end)+m);
                        end
                                
                        if normal==1 %C's spherical harmonic coefficients for the functionals with the normal field
                            if m==0
                                CElm=CEl(index((m+m_min-nminGGM):end)+m);
                            end
                        end
                            
                        if geoid==1
                            HCm=HC(index((m+m_min-nminGGM):end)+m); 
                            HSm=HS(index((m+m_min-nminGGM):end)+m);
                        end

                        Sm=S(index((m+m_min-nminGGM):end)+m);
                        %====================================================== 


                        %% Computation of the modified fnALFs
                        %======================================================
                        if volbaALFs==1 %Standard forward column method
                            if m==0
                                Pnm(:,1)=1;
                            elseif m==1                    
                                Pnm(:,1)=sqrt(3)*u.*q;  
                            elseif m>1                            
                                i=2*(2:m);
                                i1=sqrt((i+ones(size(i)))./i);
                                Pnm(:,1)=u.^m*sqrt(3)*prod(i1).*q.^m;
                            end

                            if m==nmax
                            elseif m<=(nmax-1)
                                n=m+1;
                                anm=sqrt((2*n-1)*(2*n+1)/((n-m)*(n+m)));
                                Pnm(:,2)=anm*t.*Pnm(:,1).*q;
                            end

                            if m<(nmax-1)
                                j=3;
                                for n=m+2:nmax                            
                                    anm=sqrt((2*n-1)*(2*n+1)/((n-m)*(n+m)));
                                    bnm=sqrt((2*n+1)*(n+m-1)*(n-m-1)/((n-m)*(n+m)*(2*n-3)));
                                    Pnm(:,j)=anm*t.*Pnm(:,j-1).*q-bnm*Pnm(:,j-2).*q2;
                                    j=j+1;
                                end
                            end

                        elseif volbaALFs==2 %Modified forward column method
                            if m==0
                                Pnm(:,1)=1e-280;
                            elseif m==1

                                Pnm(:,1)=sqrt(3)*q*1e-280;  
                            elseif m>1                            
                                i=2*(2:m);
                                i1=sqrt((i+ones(size(i)))./i);
                                Pnm(:,1)=sqrt(3)*prod(i1)*(q.^m)*1e-280;
                            end

                            if m==nmax
                            elseif m<=(nmax-1)
                                n=m+1;
                                anm=sqrt((2*n-1)*(2*n+1)/((n-m)*(n+m)));
                                Pnm(:,2)=anm*t.*Pnm(:,1).*q;
                            end

                            if m<(nmax-1)
                                j=3;
                                for n=m+2:nmax                            
                                    anm=sqrt((2*n-1)*(2*n+1)/((n-m)*(n+m)));
                                    bnm=sqrt((2*n+1)*(n+m-1)*(n-m-1)/((n-m)*(n+m)*(2*n-3)));
                                    Pnm(:,j)=anm*t.*Pnm(:,j-1).*q-bnm*Pnm(:,j-2).*q2;
                                    j=j+1;
                                end
                            end

                        elseif volbaALFs==3 %Extended-range arithmetic 
                            if bit==32 %32 bit version of Matlab

                                am(m+1)=rr(2*m+3);
                                for n=m+2:nmax
                                    ww=rr(2*n+1)*ri(n-m)*ri(n+m);
                                    am(n)=rr(2*n-1)*ww;
                                    bm(n)=rr(n-m-1)*rr(n+m-1)*ri(2*n-3)*ww;
                                end

                                if m~=0
                                    for i=1:length_fiG_ALFs 
                                        x=ps1(i,m);
                                        ix=ips1(i,m);

                                        if ix==0
                                            pm(m)=x;
                                        elseif ix<-1
                                            pm(m)=0;  
                                        elseif ix<0
                                            pm(m)=x*BIGI;
                                        else
                                            pm(m)=x*BIG;
                                        end

                                        if m==nmax
                                            Pnm(i,1:(nmax-m+1))=pm(m:end);
                                            continue;
                                        end

                                        y=x;
                                        iy=ix;
                                        x=(am(m+1)*t(i)*q(i))*y;
                                        ix=iy;
                                        w=abs(x);

                                        if w>=BIGS
                                            x=x*BIGI;
                                            ix=ix+1;
                                        elseif w<BIGSI
                                            x=x*BIG;
                                            ix=ix-1;
                                        end

                                        if ix==0
                                            pm(m+1)=x;
                                        elseif ix<-1  
                                            pm(m+1)=0;    
                                        elseif ix<0
                                            pm(m+1)=x*BIGI;
                                        else
                                            pm(m+1)=x*BIG;
                                        end

                                        for n=m+2:nmax 
                                            id=ix-iy;

                                            if id==0
                                                zz=(am(n)*t(i)*q(i))*x-bm(n)*q2(i)*y;
                                                iz=ix;
                                            elseif id==1
                                                zz=(am(n)*t(i)*q(i))*x-bm(n)*q2(i)*(y*BIGI);
                                                iz=ix;
                                            elseif id==-1
                                                zz=(am(n)*t(i)*q(i))*(x*BIGI)-bm(n)*q2(i)*y;
                                                iz=iy;
                                            elseif id>1
                                                zz=(am(n)*t(i)*q(i))*x;
                                                iz=ix;
                                            else
                                                zz=-bm(n)*q2(i)*y;
                                                iz=iy;
                                            end

                                            w=abs(zz);

                                            if w>=BIGS
                                                zz=zz*BIGI;
                                                iz=iz+1;
                                            elseif w<BIGSI
                                                zz=zz*BIG;
                                                iz=iz-1;
                                            end

                                            if iz==0
                                                pm(n)=zz;
                                            elseif iz<-1
                                                pm(n)=0;     
                                            elseif iz<0
                                                pm(n)=zz*BIGI;
                                            else
                                                pm(n)=zz*BIG;
                                            end

                                            y=x;
                                            iy=ix;
                                            x=zz;
                                            ix=iz;                                           
                                        end

                                        Pnm(i,1:(nmax-m+1))=pm(m:end);  
                                    end

                                elseif m==0
                                    Pnm(1:length_fiG_ALFs,1)=1;
                                    Pnm(1:length_fiG_ALFs,2)=sqrt(3)*t.*q;

                                    for i=2:nmax
                                        Pnm(1:length_fiG_ALFs,i+1)=Pnm(1:length_fiG_ALFs,i).*sqrt((2*i+1)*(2*i-1))./i.*t.*q-q2.*Pnm(1:length_fiG_ALFs,i-1).*(i-1).*sqrt(2.*i+1)./(i.*sqrt(2.*i-3));
                                    end

                                    clear rr ri am bm pm ps1 ips1 m1 m2 ...
                                        dd ix x y iy w iz zz
                                end
                            elseif bit==64 %64 bit version of Matlab

                                am(m+1)=rr(2*m+3);
                                for n=m+2:nmax
                                    ww=rr(2*n+1)*ri(n-m)*ri(n+m);
                                    am(n)=rr(2*n-1)*ww;
                                    bm(n)=rr(n-m-1)*rr(n+m-1)*ri(2*n-3)*ww;
                                end

                                if m==0 %Zonal modified fnALFs
                                    Pnm(1:length_fiG_ALFs,1)=1;
                                    Pnm(1:length_fiG_ALFs,2)=sqrt(3)*tq;
                                    for i=2:nmax
                                        Pnm(1:length_fiG_ALFs,i+1)=Pnm(1:length_fiG_ALFs,i).*sqrt((2*i+1)*(2*i-1))./i.*tq-q2.*Pnm(1:length_fiG_ALFs,i-1).*(i-1).*sqrt(2.*i+1)./(i.*sqrt(2.*i-3));
                                    end

                                    clear rr ri am bm ps1b ips1b m1 m2 dd ...
                                        ixb xb yb iyb wb izb zzb pmxb pm0b pmxBIGIb ...
                                        pmxBIGb wb wBIGSb wBIGSIb pm1xb pm10b ...
                                        pm1xBIGIb pm1xBIGb idb id0b id1b id_1b ...
                                        idv1b idm1b iz0b izm_1b izm0b izv0b tq temp1...
                                        temp2 temp3 temp4 temp5

                                elseif m~=0 %Non-zonal modified fnALFs                                    
                                    xb=ps1b(:,m);
                                    ixb=ips1b(:,m);

                                    temp5(ixb==0)=1;
                                    temp5(ixb<-1)=0;
                                    %temp5(izb>=-1 & izb<0)=BIGI;
                                    %The condition "izb>=-1 & izb<0"
                                    %is useless, as "izb" is already
                                    %initialized as "izb=BIGI".
                                    temp5(ixb>0)=BIG;
                                        
                                    Pnm(1:length_fiG_ALFs,1)=xb.*temp5;
                                    temp5=temp5.*0+BIGI; 

                                    if m<nmax
                                       yb=xb;
                                       iyb=ixb;

                                       xb=(am(m+1).*tq).*yb;
                                       ixb=iyb;
                                       wb=abs(xb);

                                       wBIGSb=wb>=BIGS;
                                       wBIGSIb=wb<BIGSI;
                                       temp3(wBIGSb)=BIGI;
                                       temp3(wBIGSIb)=BIG;
                                       temp4(wBIGSb)=1;
                                       temp4(wBIGSIb)=-1;

                                       xb=xb.*temp3;
                                       ixb=ixb+temp4;
                                       temp3=temp2;
                                       temp4=temp4.*0;

                                       temp5(ixb==0)=1;
                                       temp5(ixb<-1)=0;
                                       %temp5(izb>=-1 & izb<0)=BIGI;
                                       %The condition "izb>=-1 & izb<0"
                                       %is useless, as "izb" is already
                                       %initialized as "izb=BIGI".
                                       temp5(ixb>0)=BIG;

                                       Pnm(1:length_fiG_ALFs,2)=xb.*temp5;
                                       temp5=temp5.*0+BIGI; 

                                       for n=m+2:nmax
                                           idb=ixb-iyb;

                                           id0b=idb==0;
                                           id1b=idb==1;
                                           id_1b=idb==-1;
                                           idv1b=idb>1;

                                           temp1(id0b)=1;
                                           temp1(id1b)=1;
                                           temp2(id1b)=BIGI;
                                           temp1(id_1b)=BIGI;
                                           temp1(idv1b)=1;
                                           temp2(idv1b)=0;

                                           zzb=(am(n).*tq).*(xb.*temp1)-bm(n).*((yb.*q2).*temp2);
                                           izb=iyb;
                                           id0b_id1b_idv1b=id0b | id1b | idv1b;
                                           izb(id0b_id1b_idv1b)=ixb(id0b_id1b_idv1b);
                                           temp1=temp1.*0;
                                           temp2=temp1+1;

                                           wb=abs(zzb);

                                           wBIGSb=wb>=BIGS;
                                           wBIGSIb=wb<BIGSI;
                                           temp3(wBIGSb)=BIGI;
                                           temp3(wBIGSIb)=BIG;
                                           temp4(wBIGSb)=1;
                                           temp4(wBIGSIb)=-1;

                                           zzb=zzb.*temp3;
                                           izb=izb+temp4;
                                           temp3=temp2;
                                           temp4=temp4.*0;

                                           temp5(izb==0)=1;
                                           temp5(izb<-1)=0;
                                           %temp5(izb>=-1 & izb<0)=BIGI;
                                           %The condition "izb>=-1 & izb<0"
                                           %is useless, as "izb" is already
                                           %initialized as "izb=BIGI".
                                           temp5(izb>0)=BIG;

                                           Pnm(1:length_fiG_ALFs,n-m+1)=zzb.*temp5;
                                           temp5=temp1+BIGI;   

                                           yb=xb;
                                           iyb=ixb;
                                           xb=zzb;
                                           ixb=izb;
                                       end   
                                    end
                                end
                            end
                            
                            if symmetric_grid==true
                                if m==nmax
                                    Pnm(length_fiG_ALFs+1:end,1)=Pnm((length_fiG_ALFs-odd):-1:1,1);
                                elseif m==(nmax-1)
                                    Pnm(length_fiG_ALFs+1:end,1)=Pnm((length_fiG_ALFs-odd):-1:1,1);
                                    Pnm(length_fiG_ALFs+1:end,2)=-Pnm((length_fiG_ALFs-odd):-1:1,2);
                                else
                                    Pnm(length_fiG_ALFs+1:end,1:(nmax-m+1))=Pnm((length_fiG_ALFs-odd):-1:1,1:(nmax-m+1));
                                    Pnm(length_fiG_ALFs+1:end,2:2:(nmax-m+1))=-Pnm(length_fiG_ALFs+1:end,2:2:(nmax-m+1));
                                end
                            end
                        end
       
                        q=q_all;
                        %======================================================
                        if nmin>nminGGM
                            end_idx=min([nmin-m nmin-nminGGM]);
                            if LNOFnmin==0
                                if por==1
                                    deltaCm(1:end_idx)=0;
                                end

                                if grav==1
                                    Cm(1:end_idx)=0;
                                end

                                if geoid==1
                                    HCm(1:end_idx)=0;
                                    HSm(1:end_idx)=0;
                                end

                                Sm(1:end_idx)=0;
                            elseif LNOFnmin==1
                                Pnm(:,1:end_idx)=0;
                            end
                        end
                        %======================================================


                        %% Computation of the first-order derivatives of the modified fnALFs
                        if dALFs==1  
                            if volbaALFs==1 || volbaALFs==3
                                enm=sqrt((2*(m:nmax)+1).*((m:nmax)-m).*((m:nmax)+m)./(2*(m:nmax)-1));

                                if m==0 %Zonal modified dALFs
                                    dPnm(:,1)=0.*u;
                                    dPnm(:,2)=sqrt(3)*u.*q;
                                    dPnm(:,3:end)=-bsxfun(@times,2:nmax,tu).*Pnm(:,3:end)+(bsxfun(@times,enm(3:end),qu)).*Pnm(:,2:(end-1));
                                elseif m==nmax %Sectorial modified dALFs
                                    dPnm(:,1)=-m*(tu).*Pnm(:,1);
                                else
                                    dPnm(:,1)=-m*(tu).*Pnm(:,1); %Tesseral modified dALFs
                                    dPnm(:,2:(nmax-m+1))=-bsxfun(@times,(m+1):nmax,tu).*Pnm(:,2:(nmax-m+1))+bsxfun(@times,enm(2:end),qu).*Pnm(:,1:(nmax-m)); %Teseralne ALFs
                                end

                            elseif volbaALFs==2
                                enm=sqrt((2*(m:nmax)+1).*((m:nmax)-m).*((m:nmax)+m)./(2*(m:nmax)-1));

                                if m==0 %Zonal modified dALFs
                                    dPnm(:,1)=0.*u*1e-280;
                                    dPnm(:,2)=sqrt(3)*u.*q*1e-280;
                                    dPnm(:,3:end)=-bsxfun(@times,2:nmax,tu).*Pnm(:,3:end)+(bsxfun(@times,enm(3:end),qu)).*Pnm(:,2:(end-1));
                                elseif m==nmax
                                    dPnm(:,1)=-m*(tu).*Pnm(:,1); %Sectorial modified dALFs
                                else
                                    dPnm(:,1)=-m*(tu).*Pnm(:,1); %Tesseral modified dALFs
                                    dPnm(:,2:(nmax-m+1))=-bsxfun(@times,(m+1):nmax,tu).*Pnm(:,2:(nmax-m+1))+bsxfun(@times,enm(2:end),qu).*Pnm(:,1:(nmax-m)); %Teseralne dALFs
                                end
                            end

                            %Treatment of the dALFs singularity
                            dPnm(singdPnm,:)=0;

                            if ddALFs==1 %If the second-order derivatives of the modified fnALFs are to be computed
                                
                                if m==0 %Zonal modified ddALFs
                                    ddPnm=bsxfun(@times,tu,dPnm)-bsxfun(@times,(0:nmax).*((0:nmax)+1),Pnm);
                                else
                                    ddPnm(:,1:end-m)=bsxfun(@times,tu,dPnm(:,1:end-m))+bsxfun(@times,m^2./u.^2,Pnm(:,1:end-m))-bsxfun(@times,(m:nmax).*((m:nmax)+1),Pnm(:,1:end-m));
                                end
                                                                
                                %Treatment of the ddALFs singularity
                                ddPnm(singddPnm,:)=0;
                            end                                                                               
                        end
                            
                        m_max=max([m nminGGM]);
                        
                        %% Loop for 1:NF (number of computing functionals)
                        for i=1:pocetpar                   
                            if volbapar(i)==1       
                            elseif volbapar(i)==2 %Deflection of the vertical eta                       
                                
                                Lm=m*Pnm(:,m_min:(nmax-m+1));

                                if volbaALFs==1 || volbaALFs==3 %Computation using the standard 
                                    %forward column method or extended-range arithmetic
                                    if m==nmax
                                        Aeta=zeros(length_fiG,nmax+1);
                                        Beta=zeros(length_fiG,nmax+1);
                                    end

                                    Aeta(:,m+1)=Lm*deltaCm;
                                    Beta(:,m+1)=Lm*Sm;                        
                                elseif volbaALFs==2 %Computation using the modified forward column method combined with Horner's scheme
                                    if m==nmax
                                        eta=zeros(length_fiG,length_lambda);
                                    end

                                    eta=bsxfun(@times,eta,u)+(-Lm*deltaCm*sin(m*lambda')+Lm*Sm*cos(m*lambda'));
                                end

                            elseif volbapar(i)==3 %Deflection of the vertical xi

                                if m==nmax
                                    Aksi=zeros(length_fiG,nmax+1);
                                    Bksi=zeros(length_fiG,nmax+1);
                                end
                                
                                dLm=dPnm(:,m_min:(nmax-m+1));

                                if volbaALFs==1 || volbaALFs==3
                                    if m==nmax
                                        Aksi=zeros(length_fiG,nmax+1);
                                        Bksi=zeros(length_fiG,nmax+1);
                                    end

                                    Aksi(:,m+1)=dLm*deltaCm;
                                    Bksi(:,m+1)=dLm*Sm; 
                                elseif volbaALFs==2
                                    if m==nmax
                                        ksi=zeros(length_fiG,length_lambda);
                                    end

                                    ksi=bsxfun(@times,ksi,u)+(dLm*deltaCm*cos(m*lambda')+dLm*Sm*sin(m*lambda'));
                                end

                            elseif volbapar(i)==4 %Deflection of the vertical Theta
                               
                                Lm=m*Pnm(:,m_min:(nmax-m+1));
                                dLm=dPnm(:,m_min:(nmax-m+1));

                                if volbaALFs==1 || volbaALFs==3
                                    if m==nmax
                                        ATeta=zeros(length_fiG,nmax+1);
                                        BTeta=zeros(length_fiG,nmax+1);
                                        ATksi=ATeta;
                                        BTksi=BTeta;
                                    end

                                    ATeta(:,m+1)=Lm*deltaCm;
                                    BTeta(:,m+1)=Lm*Sm;
                                    ATksi(:,m+1)=dLm*deltaCm;
                                    BTksi(:,m+1)=dLm*Sm; 
                                elseif volbaALFs==2
                                    if m==nmax
                                        Teta=zeros(length_fiG,length_lambda);
                                        Tksi=zeros(length_fiG,length_lambda);
                                    end

                                    Teta=bsxfun(@times,Teta,u)+(-Lm*deltaCm*sin(m*lambda')+Lm*Sm*cos(m*lambda'));
                                    Tksi=bsxfun(@times,Tksi,u)+(dLm*deltaCm*cos(m*lambda')+dLm*Sm*sin(m*lambda'));
                                end

                            elseif volbapar(i)==5 %Disturbing potential
                                
                                Lm=Pnm(:,m_min:(nmax-m+1));
                                
                                if volbaALFs==1 || volbaALFs==3
                                    if m==nmax
                                        AT=zeros(length_fiG,nmax+1);
                                        BT=zeros(length_fiG,nmax+1);
                                    end

                                    AT(:,m+1)=Lm*deltaCm;
                                    BT(:,m+1)=Lm*Sm;                                    
                                elseif volbaALFs==2
                                    if m==nmax
                                        T=zeros(length_fiG,length_lambda);
                                    end

                                    T=bsxfun(@times,T,u)+(Lm*deltaCm*cos(m*lambda')+Lm*Sm*sin(m*lambda'));
                                end
                                
                            elseif volbapar(i)==6 %Disturbing tensor Trr_Tpp_Tll
                                
                                if m==nmax                               
                                    ampl_Trr=((0:nmax)+1).*((0:nmax)+2);
                                end

                                Lm=bsxfun(@times,Pnm(:,m_min:(nmax-m+1)),ampl_Trr((m_max+1):end));
                                Lmff=ddPnm(:,m_min:(nmax-m+1));
                                Lmll=m^2*Pnm(:,m_min:(nmax-m+1));

                                if volbaALFs==1 || volbaALFs==3
                                    if m==nmax
                                        ATrr=zeros(length_fiG,nmax+1);
                                        BTrr=zeros(length_fiG,nmax+1);
                                        ATff=ATrr;
                                        BTff=ATrr;
                                        ATll=ATrr;
                                        BTll=ATrr;
                                    end

                                    ATrr(:,m+1)=Lm*deltaCm;
                                    BTrr(:,m+1)=Lm*Sm;  
                                    
                                    ATff(:,m+1)=Lmff*deltaCm;
                                    BTff(:,m+1)=Lmff*Sm;
                                    
                                    ATll(:,m+1)=Lmll*deltaCm;
                                    BTll(:,m+1)=Lmll*Sm;
                                elseif volbaALFs==2
                                    if m==nmax
                                        Trr=zeros(length_fiG,length_lambda);
                                        Tff=Trr;
                                        Tll=Trr;
                                    end

                                    Trr=bsxfun(@times,Trr,u)+(Lm*deltaCm*cos(m*lambda')+Lm*Sm*sin(m*lambda'));
                                    Tff=bsxfun(@times,Tff,u)+(Lmff*deltaCm*cos(m*lambda')+Lmff*Sm*sin(m*lambda'));
                                    Tll=bsxfun(@times,Tll,u)+(Lmll*deltaCm*cos(m*lambda')+Lmll*Sm*sin(m*lambda'));
                                end
                                
                            elseif volbapar(i)==7 %Disturbing tensor Trp_Trl_Tpl
                                
                                if m==nmax                               
                                    ampl_Tfl=(0:nmax)+1;
                                end
                                
                                Lmrf=bsxfun(@times,dPnm(:,m_min:(nmax-m+1)),ampl_Tfl((m_max+1):end));
                                Lmrl=bsxfun(@times,m*Pnm(:,m_min:(nmax-m+1)),ampl_Tfl((m_max+1):end));
                                Lmfl=m*dPnm(:,m_min:(nmax-m+1));

                                if volbaALFs==1 || volbaALFs==3
                                    if m==nmax
                                        ATrf=zeros(length_fiG,nmax+1);
                                        BTrf=zeros(length_fiG,nmax+1);
                                        ATrl=ATrf;
                                        BTrl=ATrf;
                                        ATfl=ATrf;
                                        BTfl=ATrf;
                                    end

                                    ATrf(:,m+1)=Lmrf*deltaCm;
                                    BTrf(:,m+1)=Lmrf*Sm;  
                                    
                                    ATrl(:,m+1)=Lmrl*deltaCm;
                                    BTrl(:,m+1)=Lmrl*Sm;
                                    
                                    ATfl(:,m+1)=Lmfl*deltaCm;
                                    BTfl(:,m+1)=Lmfl*Sm;
                                elseif volbaALFs==2
                                    if m==nmax
                                        Trf=zeros(length_fiG,length_lambda);
                                        Trl=Trf;
                                        Tfl=Trf;
                                    end

                                    Trf=bsxfun(@times,Trf,u)+(Lmrf*deltaCm*cos(m*lambda')+Lmrf*Sm*sin(m*lambda'));
                                    Trl=bsxfun(@times,Trl,u)+(Lmrl*Sm*cos(m*lambda')-Lmrl*deltaCm*sin(m*lambda'));
                                    Tfl=bsxfun(@times,Tfl,u)+(Lmfl*Sm*cos(m*lambda')-Lmfl*deltaCm*sin(m*lambda'));
                                end
                                
                            elseif volbapar(i)==8 %Disturbing tensor Txx_Tyy_Tzz

                                if m==nmax                               
                                    ampl_Tzz=((0:nmax)+1).*((0:nmax)+2);
                                end

                                n_temp=m_max:nmax;
                                
                                if m<2
                                    Lm=bsxfun(@times,Pnm(:,m_min:(nmax-m+1)),ampl_Tzz((m_max+1):end));
                                       
                                    %Txx
                                    bnm=(n_temp+m+1).*(n_temp+m+2)./2./(m+1);
                                    LmTxx1=bsxfun(@times,Pnm(:,m_min:(nmax-m+1)),bnm-(n_temp+1).*(n_temp+2));
                                    LmTyy1=bsxfun(@times,Pnm(:,m_min:(nmax-m+1)),bnm);
                                    
                                    if m==0
                                        anm=sqrt(2)/4.*sqrt(n_temp.^2-(m+2-1).^2).*sqrt(n_temp+m+2).*sqrt(n_temp-(m+2)+2);
                                    else
                                        anm=1./4.*sqrt(n_temp.^2-(m+2-1).^2).*sqrt(n_temp+m+2).*sqrt(n_temp-(m+2)+2);
                                    end
                                    
                                    LmTxx2=bsxfun(@times,Pnm(:,m_min:(nmax-m+1)),anm);
                                    LmTxx3=zeros(length_fi,nmax-1); %Coefficients cnm are equal to zeros, if m==0 a m==1
                                else
                                    Lm=bsxfun(@times,Pnm(:,m_min:(nmax-m+1)),ampl_Tzz((m_max+1):end));
                                    
                                    %Txx
                                    bnm=(n_temp.^2+m^2+3.*n_temp+2)./2;
                                    LmTxx1=bsxfun(@times,Pnm(:,m_min:(nmax-m+1)),bnm-(n_temp+1).*(n_temp+2));
                                    LmTyy1=bsxfun(@times,Pnm(:,m_min:(nmax-m+1)),bnm);

                                    anm=1/4.*sqrt(n_temp.^2-(m+2-1).^2).*sqrt(n_temp+m+2).*sqrt(n_temp-(m+2)+2);
                                    LmTxx2=bsxfun(@times,Pnm(:,m_min:(nmax-m+1)),anm);
                                       
                                    if m==2
                                        cnm=sqrt(2)/4.*sqrt(n_temp.^2-(m-2+1).^2).*sqrt(n_temp-(m-2)).*sqrt(n_temp+m-2+2);
                                    else
                                        cnm=1./4.*sqrt(n_temp.^2-(m-2+1).^2).*sqrt(n_temp-(m-2)).*sqrt(n_temp+m-2+2);
                                    end
                                    
                                    LmTxx3=bsxfun(@times,Pnm(:,m_min:(nmax-m+1)),cnm);                                    
                                end
                   
                                if m==nmax
                                    ATzz=zeros(length_fiG,nmax+1);
                                    BTzz=zeros(length_fiG,nmax+1);
                                    ATxx1=ATzz;
                                    BTxx1=ATzz;
                                    ATxx2=ATzz;
                                    BTxx2=ATzz;
                                    ATxx3=ATzz;
                                    BTxx3=BTzz;
                                    ATyy1=ATzz;
                                    BTyy1=BTzz;
                                end

                                ATzz(:,m+1)=Lm*deltaCm;
                                BTzz(:,m+1)=Lm*Sm;  
                          
                                %Txx
                                ATxx1(:,m+1)=LmTxx1*deltaCm;
                                BTxx1(:,m+1)=LmTxx1*Sm;
 
                                if m==0 && nmax<=2
                                    ATxx2(:,m+1)=LmTxx2*deltaC(index((m+m_min-nminGGM):end)+m+2);
                                    BTxx2(:,m+1)=LmTxx2*S(index((m+m_min-nminGGM):end)+m+2);

                                    ATxx3(:,m+1)=0;
                                    BTxx3(:,m+1)=0;
                                elseif m<2 && nmax<=2
                                    ATxx2(:,m+1)=0;
                                    BTxx2(:,m+1)=0;
                                    
                                    ATxx3(:,m+1)=0;
                                    BTxx3(:,m+1)=0;
                                elseif m<2
                                    ATxx2(:,m+1)=LmTxx2*deltaC(index((m+m_min-nminGGM):end)+m+2);
                                    BTxx2(:,m+1)=LmTxx2*S(index((m+m_min-nminGGM):end)+m+2);

                                    ATxx3(:,m+1)=0;
                                    BTxx3(:,m+1)=0;
                                elseif m>nmax-2
                                    ATxx2(:,m+1)=0;
                                    BTxx2(:,m+1)=0;

                                    ATxx3(:,m+1)=LmTxx3*deltaC(index((m+m_min-nminGGM):end)+m-2);
                                    BTxx3(:,m+1)=LmTxx3*S(index((m+m_min-nminGGM):end)+m-2);
                                else
                                    ATxx2(:,m+1)=LmTxx2*deltaC(index((m+m_min-nminGGM):end)+m+2);
                                    BTxx2(:,m+1)=LmTxx2*S(index((m+m_min-nminGGM):end)+m+2);

                                    ATxx3(:,m+1)=LmTxx3*deltaC(index((m+m_min-nminGGM):end)+m-2);
                                    BTxx3(:,m+1)=LmTxx3*S(index((m+m_min-nminGGM):end)+m-2);
                                end

                                %Tyy
                                ATyy1(:,m+1)=LmTyy1*deltaCm;
                                BTyy1(:,m+1)=LmTyy1*Sm;
                                    
                                %Modified forward column method cannot be
                                %applied 
                                
                            elseif volbapar(i)==9 %Disturbing tensor Txy_Txz_Tyz
                                
                                n_temp=m:nmax;
                                
                                if m<2                                       
                                    if m==0
                                        betanm=(n_temp+2)./2.*sqrt(1+ones(1,nmax+1)).*sqrt(n_temp+m+1).*sqrt(n_temp-(m+1)+1);
                                        gamanm=zeros(1,nmax+1);
                                        
                                        %Txy
                                        dnm=1/4.*sqrt((2*n_temp+1)./(2*n_temp-1)).*sqrt(2).*sqrt(n_temp.^2-(m+2-1).^2).*sqrt(n_temp+m+2).*sqrt(n_temp+m+2-2);                                      
                                        LmTxy1=bsxfun(@times,[zeros(length_fiG,1) Pnm(:,m_min:(nmax-m))],dnm);
                                        LmTxy2=zeros(length_fi,nmax+1);
                                        LmTxy3=zeros(length_fi,nmax+1);
                                        
                                        %Tyz
                                        minm=(n_temp+2)./2.*sqrt((2*n_temp+1)./(2*n_temp-1)).*sqrt(2).*sqrt(n_temp+m+1).*sqrt(n_temp+m+1-1);
                                        LmTyz1=bsxfun(@times,[zeros(length_fiG,1) Pnm(:,m_min:(nmax-m))],minm);                                       
                                        LmTyz2=zeros(length_fi,nmax+1);
                                    else
                                        betanm=(n_temp+2)./2.*sqrt(n_temp+m+1).*sqrt(n_temp-(m+1)+1);
                                        gamanm=-(n_temp+2).*sqrt(n_temp.*(n_temp+1)./2);
                                        
                                        %Txy
                                        dnm=1/4.*sqrt((2*n_temp+1)./(2*n_temp-1)).*sqrt(n_temp.^2-(m+2-1).^2).*sqrt(n_temp+m+2).*sqrt(n_temp+m+2-2);                                        
                                        LmTxy1=bsxfun(@times,[zeros(length_fiG,1) Pnm(:,m_min:(nmax-m))],dnm);
                                        
                                        gnm=-1/4*sqrt((2*n_temp+1)./(2*n_temp-1)).*sqrt(n_temp+1).*sqrt(n_temp-1).*(n_temp+2);
                                        LmTxy2=bsxfun(@times,[zeros(length_fiG,1) Pnm(:,m_min:(nmax-m))],gnm);                                     
                                        LmTxy3=zeros(length_fi,nmax);
                                        
                                        %Tyz
                                        minm=(n_temp+2)./2.*sqrt((2*n_temp+1)./(2*n_temp-1)).*sqrt(n_temp+m+1).*sqrt(n_temp+m+1-1);
                                        LmTyz1=bsxfun(@times,[zeros(length_fiG,1) Pnm(:,m_min:(nmax-m))],minm);                                        
                                        LmTyz2=zeros(length_fi,nmax);  
                                    end
                                    
                                    %Txz
                                    LmTxz1=bsxfun(@times,Pnm(:,m_min:(nmax-m+1)),betanm);
                                    LmTxz2=bsxfun(@times,Pnm(:,m_min:(nmax-m+1)),gamanm); 
                                else
                                    %Txz
                                    betanm=(n_temp+2)./2.*sqrt(n_temp+m+1).*sqrt(n_temp-(m+1)+1);
                                    gamanm=-(n_temp+2)./2.*sqrt(n_temp-(m-1)).*sqrt(n_temp+m-1+1);
                                                                                                                                               
                                    LmTxz1=bsxfun(@times,Pnm(:,m_min:(nmax-m+1)),betanm);
                                    LmTxz2=bsxfun(@times,Pnm(:,m_min:(nmax-m+1)),gamanm);
                                                                     
                                    %Txy
                                    dnm=1/4.*sqrt((2*n_temp+1)./(2*n_temp-1)).*sqrt(n_temp.^2-(m+2-1).^2).*sqrt(n_temp+m+2).*sqrt(n_temp+m+2-2);
                                    LmTxy1=bsxfun(@times,[zeros(length_fiG,1) Pnm(:,m_min:(nmax-m))],dnm);
                                    
                                    gnm=-m/2*sqrt((2*n_temp+1)./(2*n_temp-1)).*sqrt(n_temp+m).*sqrt(n_temp-m);
                                    LmTxy2=bsxfun(@times,[zeros(length_fiG,1) Pnm(:,m_min:(nmax-m))],gnm);
                                    
                                    if m==2
                                        %Txy
                                        LmTxy3=zeros(length_fi,nmax-1);    
                                    elseif m==3
                                        %Txy
                                        hnm=-1/4*sqrt((2*n_temp+1)./(2*n_temp-1)).*sqrt(n_temp-3).*sqrt(n_temp-2).*sqrt(n_temp-1).*sqrt(n_temp+2);
                                        LmTxy3=bsxfun(@times,[zeros(length_fiG,1) Pnm(:,m_min:(nmax-m))],hnm);
                                    else
                                        %Txy
                                        hnm=-1/4*sqrt((2*n_temp+1)./(2*n_temp-1)).*sqrt(n_temp.^2-(m-2+1).^2).*sqrt(n_temp-(m-2)).*sqrt(n_temp-(m-2)-2);
                                        LmTxy3=bsxfun(@times,[zeros(length_fiG,1) Pnm(:,m_min:(nmax-m))],hnm);
                                    end
                                    
                                    %Tyz
                                    minm=(n_temp+2)./2.*sqrt((2*n_temp+1)./(2*n_temp-1)).*sqrt(n_temp+m+1).*sqrt(n_temp+m+1-1);
                                    LmTyz1=bsxfun(@times,[zeros(length_fiG,1) Pnm(:,m_min:(nmax-m))],minm);

                                    ninm=(n_temp+2)./2.*sqrt((2*n_temp+1)./(2*n_temp-1)).*sqrt(n_temp-(m-1)).*sqrt(n_temp-(m-1)-1);
                                    LmTyz2=bsxfun(@times,[zeros(length_fiG,1) Pnm(:,m_min:(nmax-m))],ninm);
                                end

                                if m==nmax
                                    ATxz1=zeros(length_fiG,nmax+1);
                                    BTxz1=zeros(length_fiG,nmax+1);
                                    ATxz2=ATxz1;
                                    BTxz2=ATxz1;
                                    ATxy1=ATxz1;
                                    BTxy1=ATxz1;
                                    ATxy2=ATxz1;
                                    BTxy2=ATxz1;
                                    ATxy3=ATxz1;
                                    BTxy3=ATxz1;
                                    ATyz1=ATxz1;
                                    BTyz1=ATxz1;
                                    ATyz2=ATxz1;
                                    BTyz2=ATxz1;
                                end

                                ATxy2(:,m+1)=(LmTxy2*deltaCm).*q;
                                BTxy2(:,m+1)=(LmTxy2*Sm).*q;
                                
                                if m==0 && nmax<=2
                                    ATxz1(:,m+1)=LmTxz1*deltaC(index((m+m_min-nminGGM):end)+m+1);
                                    BTxz1(:,m+1)=LmTxz1*S(index((m+m_min-nminGGM):end)+m+1);
                                    
                                    ATxz2(:,m+1)=0;
                                    BTxz2(:,m+1)=0;
                                    
                                    ATyz2(:,m+1)=0;
                                    BTyz2(:,m+1)=0;
                                    
                                    ATxy1(:,m+1)=(LmTxy1*deltaC(index((m+m_min-nminGGM):end)+m+2)).*q;
                                    BTxy1(:,m+1)=(LmTxy1*S(index((m+m_min-nminGGM):end)+m+2)).*q;
                                    
                                    ATxy3(:,m+1)=0;
                                    BTxy3(:,m+1)=0;
                                    
                                    ATyz1(:,m+1)=(LmTyz1*deltaC(index((m+m_min-nminGGM):end)+m+1)).*q; 
                                    BTyz1(:,m+1)=(LmTyz1*S(index((m+m_min-nminGGM):end)+m+1)).*q; 
                                    
                                elseif m<2 && nmax<=2
                                    
                                    ATxz1(:,m+1)=LmTxz1*deltaC(index((m+m_min-nminGGM):end)+m+1);
                                    BTxz1(:,m+1)=LmTxz1*S(index((m+m_min-nminGGM):end)+m+1);
                                    
                                    ATxz2(:,m+1)=LmTxz2*deltaC(index((m+m_min-nminGGM):end)+m-1);
                                    BTxz2(:,m+1)=LmTxz2*S(index((m+m_min-nminGGM):end)+m-1);
                                    
                                    ATyz2(:,m+1)=LmTyz2*deltaC(index((m+m_min-nminGGM):end)+m-1).*q;
                                    BTyz2(:,m+1)=LmTyz2*S(index((m+m_min-nminGGM):end)+m-1).*q;
                                    
                                    ATxy1(:,m+1)=0;
                                    BTxy1(:,m+1)=0;
                                    
                                    ATxy3(:,m+1)=0;
                                    BTxy3(:,m+1)=0;
                                    
                                    ATyz1(:,m+1)=(LmTyz1*deltaC(index((m+m_min-nminGGM):end)+m+1)).*q; 
                                    BTyz1(:,m+1)=(LmTyz1*S(index((m+m_min-nminGGM):end)+m+1)).*q; 

                                elseif m<2
                                    ATxz1(:,m+1)=LmTxz1*deltaC(index((m+m_min-nminGGM):end)+m+1);
                                    BTxz1(:,m+1)=LmTxz1*S(index((m+m_min-nminGGM):end)+m+1);

                                    if m==1
                                        %Txz
                                        ATxz2(:,m+1)=LmTxz2*deltaC(index((m+m_min-nminGGM):end)+m-1);
                                        BTxz2(:,m+1)=LmTxz2*S(index((m+m_min-nminGGM):end)+m-1);
                                                
                                        %Tyz
                                        ATyz2(:,m+1)=LmTyz2*deltaC(index((m+m_min-nminGGM):end)+m-1).*q;
                                        BTyz2(:,m+1)=LmTyz2*S(index((m+m_min-nminGGM):end)+m-1).*q;
                                    else
                                        %Txz
                                        ATxz2(:,m+1)=0;
                                        BTxz2(:,m+1)=0;
                                                
                                        %Tyz
                                        ATyz2(:,m+1)=0;
                                        BTyz2(:,m+1)=0;
                                    end
                                        
                                    %Txy
                                    ATxy1(:,m+1)=(LmTxy1*deltaC(index((m+m_min-nminGGM):end)+m+2)).*q;  
                                    BTxy1(:,m+1)=(LmTxy1*S(index((m+m_min-nminGGM):end)+m+2)).*q; 

                                    ATxy3(:,m+1)=0;
                                    BTxy3(:,m+1)=0;
                                            
                                    %Tyz
                                    ATyz1(:,m+1)=(LmTyz1*deltaC(index((m+m_min-nminGGM):end)+m+1)).*q; 
                                    BTyz1(:,m+1)=(LmTyz1*S(index((m+m_min-nminGGM):end)+m+1)).*q; 
                                elseif m>nmax-2
                                            
                                    if m==nmax
                                        %Txz
                                        ATxz1(:,m+1)=0;
                                        BTxz1(:,m+1)=0;
                                                
                                        %Tyz
                                        ATyz1(:,m+1)=0;
                                        BTyz1(:,m+1)=0;
                                    else
                                        %Txz
                                        ATxz1(:,m+1)=LmTxz1*deltaC(index((m+m_min-nminGGM):end)+m+1);
                                        BTxz1(:,m+1)=LmTxz1*S(index((m+m_min-nminGGM):end)+m+1);
                                                
                                        %Tyz
                                        ATyz1(:,m+1)=LmTyz1*deltaC(index((m+m_min-nminGGM):end)+m+1).*q;
                                        BTyz1(:,m+1)=LmTyz1*S(index((m+m_min-nminGGM):end)+m+1).*q;
                                    end
                                           
                                    ATxz2(:,m+1)=LmTxz2*deltaC(index((m+m_min-nminGGM):end)+m-1);
                                    BTxz2(:,m+1)=LmTxz2*S(index((m+m_min-nminGGM):end)+m-1);

                                    %Txy
                                    ATxy1(:,m+1)=0;
                                    BTxy1(:,m+1)=0;
                                            
                                    ATxy3(:,m+1)=(LmTxy3*deltaC(index((m+m_min-nminGGM):end)+m-2)).*q;
                                    BTxy3(:,m+1)=(LmTxy3*S(index((m+m_min-nminGGM):end)+m-2)).*q;  
                                            
                                    %Tyz
                                    ATyz2(:,m+1)=LmTyz2*deltaC(index((m+m_min-nminGGM):end)+m-1).*q;
                                    BTyz2(:,m+1)=LmTyz2*S(index((m+m_min-nminGGM):end)+m-1).*q;
                                else
                                    ATxz1(:,m+1)=LmTxz1*deltaC(index((m+m_min-nminGGM):end)+m+1);
                                    BTxz1(:,m+1)=LmTxz1*S(index((m+m_min-nminGGM):end)+m+1);

                                    ATxz2(:,m+1)=LmTxz2*deltaC(index((m+m_min-nminGGM):end)+m-1);
                                    BTxz2(:,m+1)=LmTxz2*S(index((m+m_min-nminGGM):end)+m-1);

                                    %Txy
                                    ATxy1(:,m+1)=(LmTxy1*deltaC(index((m+m_min-nminGGM):end)+m+2)).*q;
                                    BTxy1(:,m+1)=(LmTxy1*S(index((m+m_min-nminGGM):end)+m+2)).*q;
                                            
                                    ATxy3(:,m+1)=(LmTxy3*deltaC(index((m+m_min-nminGGM):end)+m-2)).*q;
                                    BTxy3(:,m+1)=(LmTxy3*S(index((m+m_min-nminGGM):end)+m-2)).*q;
                                            
                                    %Tyz
                                    ATyz1(:,m+1)=LmTyz1*deltaC(index((m+m_min-nminGGM):end)+m+1).*q;
                                    BTyz1(:,m+1)=LmTyz1*S(index((m+m_min-nminGGM):end)+m+1).*q;
                                           
                                    ATyz2(:,m+1)=LmTyz2*deltaC(index((m+m_min-nminGGM):end)+m-1).*q;
                                    BTyz2(:,m+1)=LmTyz2*S(index((m+m_min-nminGGM):end)+m-1).*q;
                                end
                  
                                %Modified forward column method cannot be
                                %applied
                                
                            elseif volbapar(i)==10 %Geoid undulation

                                if m==nmax
                                    amplH=zeros(length_fiG,nmax+1); %Damping factor
                                    for n=0:nmax
                                        amplH(:,n+1)=1./((R./r).^n);

                                        % When computing H, there is no
                                        % dumping factor (R./r).^n,
                                        % therefore the matrix Pnm has to be
                                        % devided by 1./((R./r).^n), since
                                        % Pnm is the matrix of the MODIFIED
                                        % fnALFS
                                    end
                                end
                                        
                                Lm=Pnm(:,m_min:(nmax-m+1));
                                LmH=Pnm(:,m_min:(nmax-m+1)).*amplH(:,(m_max+1):end);
                                
                                if volbaALFs==1 || volbaALFs==3
                                    if m==nmax
                                        AN1c=zeros(length_fiG,nmax+1);
                                        BN1c=AN1c;
                                        AH=AN1c;
                                        BH=AN1c;
                                    end  
                                    
                                    AN1c(:,m+1)=Lm*deltaCm;
                                    BN1c(:,m+1)=Lm*Sm;  
                                    AH(:,m+1)=LmH*HCm;
                                    BH(:,m+1)=LmH*HSm;
                                elseif volbaALFs==2
                                    if m==nmax
                                        N1c=zeros(length_fiG,length_lambda);
                                        H=N1c;
                                    end

                                    N1c=bsxfun(@times,N1c,u)+(Lm*deltaCm*cos(m*lambda')+Lm*Sm*sin(m*lambda'));
                                    H=bsxfun(@times,H,u)+(LmH*HCm*cos(m*lambda')+LmH*HSm*sin(m*lambda'));
                                end
                                
                            elseif volbapar(i)==11 %Gravitational potential

                                Lm=Pnm(:,m_min:(nmax-m+1));

                                if volbaALFs==1 || volbaALFs==3
                                    if m==nmax
                                        AV=zeros(length_fiG,nmax+1);
                                        BV=zeros(length_fiG,nmax+1);
                                    end

                                    AV(:,m+1)=Lm*Cm;
                                    BV(:,m+1)=Lm*Sm;  
                                elseif volbaALFs==2
                                    if m==nmax
                                        V=zeros(length_fiG,length_lambda);
                                    end

                                    V=bsxfun(@times,V,u)+(Lm*Cm*cos(m*lambda')+Lm*Sm*sin(m*lambda'));
                                end

                            elseif volbapar(i)==12 %Gravitational tensor Vrr_Vpp_Vll
                                
                                if m==nmax                               
                                    ampl_Vrr=((0:nmax)+1).*((0:nmax)+2);
                                end

                                Lm=bsxfun(@times,Pnm(:,m_min:(nmax-m+1)),ampl_Vrr((m_max+1):end));
                                Lmff=ddPnm(:,m_min:(nmax-m+1));
                                Lmll=m^2*Pnm(:,m_min:(nmax-m+1));

                                if volbaALFs==1 || volbaALFs==3
                                    if m==nmax
                                        AVrr=zeros(length_fiG,nmax+1);
                                        BVrr=zeros(length_fiG,nmax+1);
                                        AVff=AVrr;
                                        BVff=AVrr;
                                        AVll=AVrr;
                                        BVll=AVrr;
                                    end

                                    AVrr(:,m+1)=Lm*Cm;
                                    BVrr(:,m+1)=Lm*Sm;  
                                    
                                    AVff(:,m+1)=Lmff*Cm;
                                    BVff(:,m+1)=Lmff*Sm;
                                    
                                    AVll(:,m+1)=Lmll*Cm;
                                    BVll(:,m+1)=Lmll*Sm;
                                elseif volbaALFs==2
                                    if m==nmax
                                        Vrr=zeros(length_fiG,length_lambda);
                                        Vff=Vrr;
                                        Vll=Vrr;
                                    end

                                    Vrr=bsxfun(@times,Vrr,u)+(Lm*Cm*cos(m*lambda')+Lm*Sm*sin(m*lambda'));
                                    Vff=bsxfun(@times,Vff,u)+(Lmff*Cm*cos(m*lambda')+Lmff*Sm*sin(m*lambda'));
                                    Vll=bsxfun(@times,Vll,u)+(Lmll*Cm*cos(m*lambda')+Lmll*Sm*sin(m*lambda'));
                                end
                                
                            elseif volbapar(i)==13 %Gravitational tensor Vrp_Vrl_Vpl
                                
                                if m==nmax                               
                                    ampl_Vfl=((0:nmax)+1);
                                end
                                
                                Lmrf=bsxfun(@times,dPnm(:,m_min:(nmax-m+1)),ampl_Vfl((m_max+1):end));
                                Lmrl=bsxfun(@times,m*Pnm(:,m_min:(nmax-m+1)),ampl_Vfl((m_max+1):end));
                                Lmfl=m*dPnm(:,m_min:(nmax-m+1));

                                if volbaALFs==1 || volbaALFs==3
                                    if m==nmax
                                        AVrf=zeros(length_fiG,nmax+1);
                                        BVrf=zeros(length_fiG,nmax+1);
                                        AVrl=AVrf;
                                        BVrl=AVrf;
                                        AVfl=AVrf;
                                        BVfl=AVrf;
                                    end

                                    AVrf(:,m+1)=Lmrf*Cm;
                                    BVrf(:,m+1)=Lmrf*Sm;  
                                    
                                    AVrl(:,m+1)=Lmrl*Cm;
                                    BVrl(:,m+1)=Lmrl*Sm;
                                    
                                    AVfl(:,m+1)=Lmfl*Cm;
                                    BVfl(:,m+1)=Lmfl*Sm;
                                elseif volbaALFs==2
                                    if m==nmax
                                        Vrf=zeros(length_fiG,length_lambda);
                                        Vrl=Vrf;
                                        Vfl=Vrf;
                                    end

                                    Vrf=bsxfun(@times,Vrf,u)+(Lmrf*Cm*cos(m*lambda')+Lmrf*Sm*sin(m*lambda'));
                                    Vrl=bsxfun(@times,Vrl,u)+(Lmrl*Sm*cos(m*lambda')-Lmrl*Cm*sin(m*lambda'));
                                    Vfl=bsxfun(@times,Vfl,u)+(Lmfl*Sm*cos(m*lambda')-Lmfl*Cm*sin(m*lambda'));
                                end
                                
                            elseif volbapar(i)==14 %Gravitational tensor Vxx_Vyy_Vzz

                                if m==nmax                               
                                    ampl_Vzz=((0:nmax)+1).*((0:nmax)+2);
                                end

                                n_temp=m_max:nmax;
                                
                                if m<2
                                    Lm=bsxfun(@times,Pnm(:,m_min:(nmax-m+1)),ampl_Vzz((m_max+1):end));
                                       
                                    %Vxx
                                    bnm=(n_temp+m+1).*(n_temp+m+2)./2./(m+1);
                                    LmVxx1=bsxfun(@times,Pnm(:,m_min:(nmax-m+1)),bnm-(n_temp+1).*(n_temp+2));
                                    LmVyy1=bsxfun(@times,Pnm(:,m_min:(nmax-m+1)),bnm);
                                    
                                    if m==0
                                        anm=sqrt(2)/4.*sqrt(n_temp.^2-(m+2-1).^2).*sqrt(n_temp+m+2).*sqrt(n_temp-(m+2)+2);
                                    else
                                        anm=1./4.*sqrt(n_temp.^2-(m+2-1).^2).*sqrt(n_temp+m+2).*sqrt(n_temp-(m+2)+2);
                                    end
                                    
                                    LmVxx2=bsxfun(@times,Pnm(:,m_min:(nmax-m+1)),anm);
                                    LmVxx3=zeros(length_fi,nmax-1); %Coefficients cnm are equal to zero, if m==0 a m==1
                                else
                                    Lm=bsxfun(@times,Pnm(:,m_min:(nmax-m+1)),ampl_Vzz((m_max+1):end));
                                 
                                    %Vxx
                                    bnm=(n_temp.^2+m^2+3.*n_temp+2)./2;
                                    LmVxx1=bsxfun(@times,Pnm(:,m_min:(nmax-m+1)),bnm-(n_temp+1).*(n_temp+2));
                                    LmVyy1=bsxfun(@times,Pnm(:,m_min:(nmax-m+1)),bnm);

                                    anm=1/4.*sqrt(n_temp.^2-(m+2-1).^2).*sqrt(n_temp+m+2).*sqrt(n_temp-(m+2)+2);
                                    LmVxx2=bsxfun(@times,Pnm(:,m_min:(nmax-m+1)),anm);
                                       
                                    if m==2
                                        cnm=sqrt(2)/4.*sqrt(n_temp.^2-(m-2+1).^2).*sqrt(n_temp-(m-2)).*sqrt(n_temp+m-2+2);
                                    else
                                        cnm=1./4.*sqrt(n_temp.^2-(m-2+1).^2).*sqrt(n_temp-(m-2)).*sqrt(n_temp+m-2+2);
                                    end
                                    
                                    LmVxx3=bsxfun(@times,Pnm(:,m_min:(nmax-m+1)),cnm);                                    
                                end
                   
                                if m==nmax
                                    AVzz=zeros(length_fiG,nmax+1);
                                    BVzz=zeros(length_fiG,nmax+1);
                                    AVxx1=AVzz;
                                    BVxx1=AVzz;
                                    AVxx2=AVzz;
                                    BVxx2=AVzz;
                                    AVxx3=AVzz;
                                    BVxx3=BVzz;
                                    AVyy1=AVzz;
                                    BVyy1=BVzz;
                                end

                                AVzz(:,m+1)=Lm*Cm;
                                BVzz(:,m+1)=Lm*Sm;  

                                %Vxx
                                AVxx1(:,m+1)=LmVxx1*Cm;
                                BVxx1(:,m+1)=LmVxx1*Sm;
 
                                if m==0 && nmax<=2
                                    AVxx2(:,m+1)=LmVxx2*C(index((m+m_min-nminGGM):end)+m+2);
                                    BVxx2(:,m+1)=LmVxx2*S(index((m+m_min-nminGGM):end)+m+2);

                                    AVxx3(:,m+1)=0;
                                    BVxx3(:,m+1)=0;
                                elseif m<2 && nmax<=2
                                    AVxx2(:,m+1)=0;
                                    BVxx2(:,m+1)=0;
                                    
                                    AVxx3(:,m+1)=0;
                                    BVxx3(:,m+1)=0;
                                elseif m<2
                                    AVxx2(:,m+1)=LmVxx2*C(index((m+m_min-nminGGM):end)+m+2);
                                    BVxx2(:,m+1)=LmVxx2*S(index((m+m_min-nminGGM):end)+m+2);

                                    AVxx3(:,m+1)=0;
                                    BVxx3(:,m+1)=0;
                                elseif m>nmax-2
                                    AVxx2(:,m+1)=0;
                                    BVxx2(:,m+1)=0;
                                           
                                    AVxx3(:,m+1)=LmVxx3*C(index((m+m_min-nminGGM):end)+m-2);
                                    BVxx3(:,m+1)=LmVxx3*S(index((m+m_min-nminGGM):end)+m-2);
                                else
                                    AVxx2(:,m+1)=LmVxx2*C(index((m+m_min-nminGGM):end)+m+2);
                                    BVxx2(:,m+1)=LmVxx2*S(index((m+m_min-nminGGM):end)+m+2);

                                    AVxx3(:,m+1)=LmVxx3*C(index((m+m_min-nminGGM):end)+m-2);
                                    BVxx3(:,m+1)=LmVxx3*S(index((m+m_min-nminGGM):end)+m-2);
                                end

                                %Vyy
                                AVyy1(:,m+1)=LmVyy1*Cm;
                                BVyy1(:,m+1)=LmVyy1*Sm;
                                    
                                %Modified forward column method cannot be
                                %applied
                                
                            elseif volbapar(i)==15 %Gravitational tensor Vxy_Vxz_Vyz
                                
                                n_temp=m:nmax;
                                
                                if m<2                                       
                                    if m==0
                                        betanm=(n_temp+2)./2.*sqrt(1+ones(1,nmax+1)).*sqrt(n_temp+m+1).*sqrt(n_temp-(m+1)+1);
                                        gamanm=zeros(1,nmax+1);
                                        
                                        %Vxy
                                        dnm=1/4.*sqrt((2*n_temp+1)./(2*n_temp-1)).*sqrt(2).*sqrt(n_temp.^2-(m+2-1).^2).*sqrt(n_temp+m+2).*sqrt(n_temp+m+2-2);                                      
                                        LmVxy1=bsxfun(@times,[zeros(length_fiG,1) Pnm(:,m_min:(nmax-m))],dnm);
                                        LmVxy2=zeros(length_fi,nmax+1);                                       
                                        LmVxy3=zeros(length_fi,nmax+1);
                                        
                                        %Vyz
                                        minm=(n_temp+2)./2.*sqrt((2*n_temp+1)./(2*n_temp-1)).*sqrt(2).*sqrt(n_temp+m+1).*sqrt(n_temp+m+1-1);
                                        LmVyz1=bsxfun(@times,[zeros(length_fiG,1) Pnm(:,m_min:(nmax-m))],minm);
                                        LmVyz2=zeros(length_fi,nmax+1);
                                    else
                                        betanm=(n_temp+2)./2.*sqrt(n_temp+m+1).*sqrt(n_temp-(m+1)+1);
                                        gamanm=-(n_temp+2).*sqrt(n_temp.*(n_temp+1)./2);
                                        
                                        %Vxy
                                        dnm=1/4.*sqrt((2*n_temp+1)./(2*n_temp-1)).*sqrt(n_temp.^2-(m+2-1).^2).*sqrt(n_temp+m+2).*sqrt(n_temp+m+2-2);                                        
                                        LmVxy1=bsxfun(@times,[zeros(length_fiG,1) Pnm(:,m_min:(nmax-m))],dnm);
                                        
                                        gnm=-1/4*sqrt((2*n_temp+1)./(2*n_temp-1)).*sqrt(n_temp+1).*sqrt(n_temp-1).*(n_temp+2);
                                        LmVxy2=bsxfun(@times,[zeros(length_fiG,1) Pnm(:,m_min:(nmax-m))],gnm);
                                        LmVxy3=zeros(length_fi,nmax);
                                        
                                        %Vyz
                                        minm=(n_temp+2)./2.*sqrt((2*n_temp+1)./(2*n_temp-1)).*sqrt(n_temp+m+1).*sqrt(n_temp+m+1-1);
                                        LmVyz1=bsxfun(@times,[zeros(length_fiG,1) Pnm(:,m_min:(nmax-m))],minm);
                                        LmVyz2=zeros(length_fi,nmax);
                                    end
                                    
                                    %Vxz
                                    LmVxz1=bsxfun(@times,Pnm(:,m_min:(nmax-m+1)),betanm);
                                    LmVxz2=bsxfun(@times,Pnm(:,m_min:(nmax-m+1)),gamanm); 
                                else
                                    %Vxz
                                    betanm=(n_temp+2)./2.*sqrt(n_temp+m+1).*sqrt(n_temp-(m+1)+1);
                                    gamanm=-(n_temp+2)./2.*sqrt(n_temp-(m-1)).*sqrt(n_temp+m-1+1);
                                                                                                                                               
                                    LmVxz1=bsxfun(@times,Pnm(:,m_min:(nmax-m+1)),betanm);
                                    LmVxz2=bsxfun(@times,Pnm(:,m_min:(nmax-m+1)),gamanm);
                                                                     
                                    %Vxy
                                    dnm=1/4.*sqrt((2*n_temp+1)./(2*n_temp-1)).*sqrt(n_temp.^2-(m+2-1).^2).*sqrt(n_temp+m+2).*sqrt(n_temp+m+2-2);
                                    LmVxy1=bsxfun(@times,[zeros(length_fiG,1) Pnm(:,m_min:(nmax-m))],dnm);
                                    
                                    gnm=-m/2*sqrt((2*n_temp+1)./(2*n_temp-1)).*sqrt(n_temp+m).*sqrt(n_temp-m);
                                    LmVxy2=bsxfun(@times,[zeros(length_fiG,1) Pnm(:,m_min:(nmax-m))],gnm);
                                    
                                    if m==2
                                        %Vxy
                                        LmVxy3=zeros(length_fi,nmax-1);    
                                    elseif m==3
                                        %Vxy
                                        hnm=-1/4*sqrt((2*n_temp+1)./(2*n_temp-1)).*sqrt(n_temp-3).*sqrt(n_temp-2).*sqrt(n_temp-1).*sqrt(n_temp+2);
                                        LmVxy3=bsxfun(@times,[zeros(length_fiG,1) Pnm(:,m_min:(nmax-m))],hnm);
                                    else
                                        %Vxy
                                        hnm=-1/4*sqrt((2*n_temp+1)./(2*n_temp-1)).*sqrt(n_temp.^2-(m-2+1).^2).*sqrt(n_temp-(m-2)).*sqrt(n_temp-(m-2)-2);
                                        LmVxy3=bsxfun(@times,[zeros(length_fiG,1) Pnm(:,m_min:(nmax-m))],hnm);
                                    end
                                    
                                    %Vyz
                                    minm=(n_temp+2)./2.*sqrt((2*n_temp+1)./(2*n_temp-1)).*sqrt(n_temp+m+1).*sqrt(n_temp+m+1-1);
                                    LmVyz1=bsxfun(@times,[zeros(length_fiG,1) Pnm(:,m_min:(nmax-m))],minm);

                                    ninm=(n_temp+2)./2.*sqrt((2*n_temp+1)./(2*n_temp-1)).*sqrt(n_temp-(m-1)).*sqrt(n_temp-(m-1)-1);
                                    LmVyz2=bsxfun(@times,[zeros(length_fiG,1) Pnm(:,m_min:(nmax-m))],ninm);
                                end

                                if m==nmax
                                    AVxz1=zeros(length_fiG,nmax+1);
                                    BVxz1=zeros(length_fiG,nmax+1);
                                    AVxz2=AVxz1;
                                    BVxz2=AVxz1;
                                    AVxy1=AVxz1;
                                    BVxy1=AVxz1;
                                    AVxy2=AVxz1;
                                    BVxy2=AVxz1;
                                    AVxy3=AVxz1;
                                    BVxy3=AVxz1;
                                    AVyz1=AVxz1;
                                    BVyz1=AVxz1;
                                    AVyz2=AVxz1;
                                    BVyz2=AVxz1;
                                end

                                AVxy2(:,m+1)=(LmVxy2*Cm).*q;
                                BVxy2(:,m+1)=(LmVxy2*Sm).*q;
                                
                                if m==0 && nmax<=2
                                    AVxz1(:,m+1)=LmVxz1*C(index((m+m_min-nminGGM):end)+m+1);
                                    BVxz1(:,m+1)=LmVxz1*S(index((m+m_min-nminGGM):end)+m+1);
                                    
                                    AVxz2(:,m+1)=0;
                                    BVxz2(:,m+1)=0;
                                    
                                    AVyz2(:,m+1)=0;
                                    BVyz2(:,m+1)=0;
                                    
                                    AVxy1(:,m+1)=(LmVxy1*C(index((m+m_min-nminGGM):end)+m+2)).*q;
                                    BVxy1(:,m+1)=(LmVxy1*S(index((m+m_min-nminGGM):end)+m+2)).*q;
                                    
                                    AVxy3(:,m+1)=0;
                                    BVxy3(:,m+1)=0;
                                    
                                    AVyz1(:,m+1)=(LmVyz1*C(index((m+m_min-nminGGM):end)+m+1)).*q; 
                                    BVyz1(:,m+1)=(LmVyz1*S(index((m+m_min-nminGGM):end)+m+1)).*q; 
                                    
                                elseif m<2 && nmax<=2
                                    
                                    AVxz1(:,m+1)=LmVxz1*C(index((m+m_min-nminGGM):end)+m+1);
                                    BVxz1(:,m+1)=LmVxz1*S(index((m+m_min-nminGGM):end)+m+1);
                                    
                                    AVxz2(:,m+1)=LmVxz2*C(index((m+m_min-nminGGM):end)+m-1);
                                    BVxz2(:,m+1)=LmVxz2*S(index((m+m_min-nminGGM):end)+m-1);
                                    
                                    AVyz2(:,m+1)=LmVyz2*C(index((m+m_min-nminGGM):end)+m-1).*q;
                                    BVyz2(:,m+1)=LmVyz2*S(index((m+m_min-nminGGM):end)+m-1).*q;
                                    
                                    AVxy1(:,m+1)=0;
                                    BVxy1(:,m+1)=0;
                                    
                                    AVxy3(:,m+1)=0;
                                    BVxy3(:,m+1)=0;
                                    
                                    AVyz1(:,m+1)=(LmVyz1*C(index((m+m_min-nminGGM):end)+m+1)).*q; 
                                    BVyz1(:,m+1)=(LmVyz1*S(index((m+m_min-nminGGM):end)+m+1)).*q; 
                                elseif m<2
                                    AVxz1(:,m+1)=LmVxz1*C(index((m+m_min-nminGGM):end)+m+1);
                                    BVxz1(:,m+1)=LmVxz1*S(index((m+m_min-nminGGM):end)+m+1);

                                    if m==1
                                        %Vxz
                                        AVxz2(:,m+1)=LmVxz2*C(index((m+m_min-nminGGM):end)+m-1);
                                        BVxz2(:,m+1)=LmVxz2*S(index((m+m_min-nminGGM):end)+m-1);
                                                
                                        %Vyz
                                        AVyz2(:,m+1)=LmVyz2*C(index((m+m_min-nminGGM):end)+m-1).*q;
                                        BVyz2(:,m+1)=LmVyz2*S(index((m+m_min-nminGGM):end)+m-1).*q;
                                    else
                                        %Vxz
                                        AVxz2(:,m+1)=0;
                                        BVxz2(:,m+1)=0;
                                                
                                        %Vyz
                                        AVyz2(:,m+1)=0;
                                        BVyz2(:,m+1)=0;
                                    end
                                        
                                    %Vxy
                                    AVxy1(:,m+1)=(LmVxy1*C(index((m+m_min-nminGGM):end)+m+2)).*q;  
                                    BVxy1(:,m+1)=(LmVxy1*S(index((m+m_min-nminGGM):end)+m+2)).*q; 
                                          
                                    AVxy3(:,m+1)=0;
                                    BVxy3(:,m+1)=0;
                                            
                                    %Vyz
                                    AVyz1(:,m+1)=(LmVyz1*C(index((m+m_min-nminGGM):end)+m+1)).*q;  
                                    BVyz1(:,m+1)=(LmVyz1*S(index((m+m_min-nminGGM):end)+m+1)).*q; 
                                elseif m>nmax-2
                                            
                                    if m==nmax
                                        %Vxz
                                        AVxz1(:,m+1)=0;
                                        BVxz1(:,m+1)=0;
                                                
                                        %Vyz
                                        AVyz1(:,m+1)=0;
                                        BVyz1(:,m+1)=0;
                                    else
                                        %Vxz
                                        AVxz1(:,m+1)=LmVxz1*C(index((m+m_min-nminGGM):end)+m+1);
                                        BVxz1(:,m+1)=LmVxz1*S(index((m+m_min-nminGGM):end)+m+1);
                                                
                                        %Vyz
                                        AVyz1(:,m+1)=LmVyz1*C(index((m+m_min-nminGGM):end)+m+1).*q;
                                        BVyz1(:,m+1)=LmVyz1*S(index((m+m_min-nminGGM):end)+m+1).*q;
                                    end
                                           
                                    AVxz2(:,m+1)=LmVxz2*C(index((m+m_min-nminGGM):end)+m-1);
                                    BVxz2(:,m+1)=LmVxz2*S(index((m+m_min-nminGGM):end)+m-1);

                                    %Vxy
                                    AVxy1(:,m+1)=0;
                                    BVxy1(:,m+1)=0;
                                            
                                    AVxy3(:,m+1)=(LmVxy3*C(index((m+m_min-nminGGM):end)+m-2)).*q;
                                    BVxy3(:,m+1)=(LmVxy3*S(index((m+m_min-nminGGM):end)+m-2)).*q;  
                                           
                                    %Vyz
                                    AVyz2(:,m+1)=LmVyz2*C(index((m+m_min-nminGGM):end)+m-1).*q;
                                    BVyz2(:,m+1)=LmVyz2*S(index((m+m_min-nminGGM):end)+m-1).*q;
                                else
                                    AVxz1(:,m+1)=LmVxz1*C(index((m+m_min-nminGGM):end)+m+1);
                                    BVxz1(:,m+1)=LmVxz1*S(index((m+m_min-nminGGM):end)+m+1);

                                    AVxz2(:,m+1)=LmVxz2*C(index((m+m_min-nminGGM):end)+m-1);
                                    BVxz2(:,m+1)=LmVxz2*S(index((m+m_min-nminGGM):end)+m-1);

                                    %Vxy
                                    AVxy1(:,m+1)=(LmVxy1*C(index((m+m_min-nminGGM):end)+m+2)).*q;
                                    BVxy1(:,m+1)=(LmVxy1*S(index((m+m_min-nminGGM):end)+m+2)).*q;
                                           
                                    AVxy3(:,m+1)=(LmVxy3*C(index((m+m_min-nminGGM):end)+m-2)).*q;
                                    BVxy3(:,m+1)=(LmVxy3*S(index((m+m_min-nminGGM):end)+m-2)).*q;
                                            
                                    %Vyz
                                    AVyz1(:,m+1)=LmVyz1*C(index((m+m_min-nminGGM):end)+m+1).*q;
                                    BVyz1(:,m+1)=LmVyz1*S(index((m+m_min-nminGGM):end)+m+1).*q;
                                           
                                    AVyz2(:,m+1)=LmVyz2*C(index((m+m_min-nminGGM):end)+m-1).*q;
                                    BVyz2(:,m+1)=LmVyz2*S(index((m+m_min-nminGGM):end)+m-1).*q;
                                end
                                
                                %Modified forward column method cannot be
                                %applied
                                
                            elseif volbapar(i)==16 %Gravity vector gX_gY_gZ

                                if m==nmax                                
                                    ampl_Wr=((0:nmax)+1);
                                end

                                LmWr=bsxfun(@times,Pnm(:,m_min:(nmax-m+1)),ampl_Wr((m_max+1):end));
                                LmWlambda=m*Pnm(:,m_min:(nmax-m+1));
                                LmWfi=dPnm(:,m_min:(nmax-m+1));
                       

                                if volbaALFs==1 || volbaALFs==3
                                    if m==nmax
                                        AWr=zeros(length_fiG,nmax+1);
                                        BWr=zeros(length_fiG,nmax+1);
                                        AWlambda=AWr;
                                        BWlambda=BWr;
                                        AWfi=AWr;
                                        BWfi=BWr;
                                    end

                                    AWr(:,m+1)=LmWr*Cm;
                                    BWr(:,m+1)=LmWr*Sm;   
                                    AWlambda(:,m+1)=LmWlambda*Cm;
                                    BWlambda(:,m+1)=LmWlambda*Sm;
                                    AWfi(:,m+1)=LmWfi*Cm;
                                    BWfi(:,m+1)=LmWfi*Sm;
                                elseif volbaALFs==2
                                    if m==nmax
                                        Wr=zeros(length_fiG,length_lambda);
                                        Wlambda=Wr;
                                        Wfi=Wr;
                                    end

                                    Wr=bsxfun(@times,Wr,u)+(LmWr*Cm*cos(m*lambda')+LmWr*Sm*sin(m*lambda'));
                                    Wlambda=bsxfun(@times,Wlambda,u)+(-LmWlambda*Cm*sin(m*lambda')+LmWlambda*Sm*cos(m*lambda'));
                                    Wfi=bsxfun(@times,Wfi,u)+(LmWfi*Cm*cos(m*lambda')+LmWfi*Sm*sin(m*lambda'));
                                end

                            elseif volbapar(i)==17 %Gravity sa

                                if m==nmax                                
                                    ampl_g_sa=((0:nmax)+1);
                                end

                                Lm=bsxfun(@times,Pnm(:,m_min:(nmax-m+1)),ampl_g_sa((m_max+1):end));

                                if volbaALFs==1 || volbaALFs==3
                                    if m==nmax
                                        Ag_sa=zeros(length_fiG,nmax+1);
                                        Bg_sa=zeros(length_fiG,nmax+1);
                                    end

                                    Ag_sa(:,m+1)=Lm*Cm;
                                    Bg_sa(:,m+1)=Lm*Sm;                             
                                elseif volbaALFs==2
                                    if m==nmax
                                       g_sa=zeros(length_fiG,length_lambda); 
                                    end

                                    g_sa=bsxfun(@times,g_sa,u)+(Lm*Cm*cos(m*lambda')+Lm*Sm*sin(m*lambda'));
                                end

                            elseif volbapar(i)==18 %Gravity potential                           

                                Lm=Pnm(:,m_min:(nmax-m+1));

                                if volbaALFs==1 || volbaALFs==3                               
                                    if m==nmax
                                        AW=zeros(length_fiG,nmax+1);
                                        BW=zeros(length_fiG,nmax+1);
                                    end

                                    AW(:,m+1)=Lm*Cm;
                                    BW(:,m+1)=Lm*Sm;                    
                                elseif volbaALFs==2
                                    if m==nmax
                                        W=zeros(length_fiG,length_lambda);
                                    end

                                    W=bsxfun(@times,W,u)+(Lm*Cm*cos(m*lambda')+Lm*Sm*sin(m*lambda'));
                                end

                            elseif volbapar(i)==19 %Gravity anomaly sa

                                if m==nmax
                                    ampl_anomalia_sa=((0:nmax)-1);
                                end

                                Lm=bsxfun(@times,Pnm(:,m_min:(nmax-m+1)),ampl_anomalia_sa((m_max+1):end));

                                if volbaALFs==1 || volbaALFs==3
                                    if m==nmax
                                        Aanomalia_sa=zeros(length_fiG,nmax+1);
                                        Banomalia_sa=zeros(length_fiG,nmax+1);
                                    end

                                    Aanomalia_sa(:,m+1)=Lm*deltaCm;
                                    Banomalia_sa(:,m+1)=Lm*Sm;
                                elseif volbaALFs==2
                                    if m==nmax
                                        anomalia_sa=zeros(length_fiG,length_lambda);
                                    end

                                    anomalia_sa=bsxfun(@times,anomalia_sa,u)+(Lm*deltaCm*cos(m*lambda')+Lm*Sm*sin(m*lambda'));
                                end

                            elseif volbapar(i)==20 %Gravity disturbance

                                if m==nmax                                
                                    ampl_Wrpor=((0:nmax)+1);
                                end

                                LmWrpor=bsxfun(@times,Pnm(:,m_min:(nmax-m+1)),ampl_Wrpor((m_max+1):end));
                                LmWlambdapor=m*Pnm(:,m_min:(nmax-m+1));
                                LmWfipor=dPnm(:,m_min:(nmax-m+1));
                       
                                if volbaALFs==1 || volbaALFs==3
                                    if m==nmax
                                        AWrpor=zeros(length_fiG,nmax+1);
                                        BWrpor=zeros(length_fiG,nmax+1);
                                        AWlambdapor=AWrpor;
                                        BWlambdapor=BWrpor;
                                        AWfipor=AWrpor;
                                        BWfipor=BWrpor;
                                    end

                                    AWrpor(:,m+1)=LmWrpor*Cm;
                                    BWrpor(:,m+1)=LmWrpor*Sm;   
                                    AWlambdapor(:,m+1)=LmWlambdapor*Cm;
                                    BWlambdapor(:,m+1)=LmWlambdapor*Sm;
                                    AWfipor(:,m+1)=LmWfipor*Cm;
                                    BWfipor(:,m+1)=LmWfipor*Sm;

                                    if m==0
                                        AUr(:,m+1)=LmWrpor*CElm;
                                        AUfi(:,m+1)=LmWfipor*CElm;
                                    end
                                elseif volbaALFs==2
                                    if m==nmax
                                        Wrpor=zeros(length_fiG,length_lambda);
                                        Wlambdapor=Wrpor;
                                        Wfipor=Wrpor;
                                        Ur=Wrpor;
                                        Ufi=Ur;
                                    end

                                    Wrpor=bsxfun(@times,Wrpor,u)+(LmWrpor*Cm*cos(m*lambda')+LmWrpor*Sm*sin(m*lambda'));
                                    Wlambdapor=bsxfun(@times,Wlambdapor,u)+(-LmWlambdapor*Cm*sin(m*lambda')+LmWlambdapor*Sm*cos(m*lambda'));
                                    Wfipor=bsxfun(@times,Wfipor,u)+(LmWfipor*Cm*cos(m*lambda')+LmWfipor*Sm*sin(m*lambda'));

                                    if m==0                                   
                                        Ur=LmWrpor*CElm*cos(m*lambda');
                                        Ufi=LmWfipor*CElm*cos(m*lambda');
                                    end
                                end

                            elseif volbapar(i)==21 %Gravity disturbance sa

                                if m==nmax                                
                                    ampl_porucha_sa=((0:nmax)+1);
                                end

                                Lm=bsxfun(@times,Pnm(:,m_min:(nmax-m+1)),ampl_porucha_sa((m_max+1):end));

                                if volbaALFs==1 || volbaALFs==3
                                    if m==nmax
                                        Aporucha_sa=zeros(length_fiG,nmax+1);
                                        Bporucha_sa=zeros(length_fiG,nmax+1);
                                    end

                                    Aporucha_sa(:,m+1)=Lm*deltaCm;
                                    Bporucha_sa(:,m+1)=Lm*Sm;                            
                                elseif volbaALFs==2
                                    if m==nmax
                                        porucha_sa=zeros(length_fiG,length_lambda);
                                    end

                                    porucha_sa=bsxfun(@times,porucha_sa,u)+(Lm*deltaCm*cos(m*lambda')+Lm*Sm*sin(m*lambda'));
                                end

                            elseif volbapar(i)==22 %Height anomaly ell

                                Lm=Pnm(:,m_min:(nmax-m+1));

                                if volbaALFs==1 || volbaALFs==3
                                    if m==nmax
                                        AzetaEl=zeros(length_fiG,nmax+1);
                                        BzetaEl=zeros(length_fiG,nmax+1);
                                    end

                                    AzetaEl(:,m+1)=Lm*deltaCm;
                                    BzetaEl(:,m+1)=Lm*Sm;                            
                                elseif volbaALFs==2
                                    if m==nmax
                                        zetaEl=zeros(length_fiG,length_lambda);
                                    end

                                    zetaEl=bsxfun(@times,zetaEl,u)+(Lm*deltaCm*cos(m*lambda')+Lm*Sm*sin(m*lambda'));
                                end

                            elseif volbapar(i)==23 %Height anomaly
                                
                                if m==nmax
                                    ampl_zeta_H=zeros(length_fiG,nmax+1); %Damping factor
                                    for n=0:nmax
                                        ampl_zeta_H(:,n+1)=1./((R./r).^n); 
                                        % When computing H, there is no
                                        % dumping factor (R./r).^n,
                                        % therefore the matrix Pnm has to be
                                        % devided by 1./((R./r).^n), since
                                        % Pnm is the matrix of the MODIFIED
                                        % fnALFS
                                    end
                                    
                                    ampl_zeta_dg=((0:nmax)+1);
                                end

                                Lmdg=bsxfun(@times,Pnm(:,m_min:(nmax-m+1)),ampl_zeta_dg(:,(m_max+1):end));
                                Lm=Pnm(:,m_min:(nmax-m+1));
                                LmH=Pnm(:,m_min:(nmax-m+1)).*ampl_zeta_H(:,(m_max+1):end);
                                
                                if volbaALFs==1 || volbaALFs==3
                                    if m==nmax
                                        AN1c_zeta=zeros(length_fiG,nmax+1);
                                        BN1c_zeta=AN1c_zeta;
                                        Azetadg=AN1c_zeta;
                                        Bzetadg=AN1c_zeta;
                                        AH_zeta=AN1c_zeta;
                                        BH_zeta=AN1c_zeta;
                                        Azeta=AN1c_zeta;
                                        Bzeta=AN1c_zeta;
                                    end  
                                    
                                    AN1c_zeta(:,m+1)=Lm*deltaCm;
                                    BN1c_zeta(:,m+1)=Lm*Sm;  
                                    Azetadg(:,m+1)=Lmdg*deltaCm;
                                    Bzetadg(:,m+1)=Lmdg*Sm;
                                    AH_zeta(:,m+1)=LmH*HCm;
                                    BH_zeta(:,m+1)=LmH*HSm;
                                    Azeta(:,m+1)=Lm*deltaCm;
                                    Bzeta(:,m+1)=Lm*Sm;
                                elseif volbaALFs==2
                                    if m==nmax
                                        zeta_N1c=zeros(length_fiG,length_lambda);
                                        zeta_H=zeta_N1c;
                                        zeta_dg=zeta_N1c;
                                        zeta_zetaEl=zeta_N1c;
                                    end

                                    zeta_N1c=bsxfun(@times,zeta_N1c,u)+(Lm*deltaCm*cos(m*lambda')+Lm*Sm*sin(m*lambda'));
                                    zeta_H=bsxfun(@times,zeta_H,u)+(LmH*HCm*cos(m*lambda')+LmH*HSm*sin(m*lambda'));
                                    zeta_dg=bsxfun(@times,zeta_dg,u)+(Lmdg*deltaCm*cos(m*lambda')+Lmdg*Sm*sin(m*lambda'));
                                    zeta_zetaEl=bsxfun(@times,zeta_zetaEl,u)+(Lm*deltaCm*cos(m*lambda')+Lm*Sm*sin(m*lambda'));
                                end
                                
                                if m==0
                                    clear ampl_zeta_H ampl_zeta_dg
                                end
                                                                
                            elseif volbapar(i)==24 %Second radial derivative of disturbing potential

                                if m==nmax                               
                                    ampl_T_rr=((0:nmax)+1).*((0:nmax)+2);
                                end
                                
                                Lm=bsxfun(@times,Pnm(:,m_min:(nmax-m+1)),ampl_T_rr(:,(m_max+1):end));
                                
                                if volbaALFs==1 || volbaALFs==3
                                    if m==nmax
                                        AT_rr=zeros(length_fiG,nmax+1);
                                        BT_rr=zeros(length_fiG,nmax+1);
                                    end

                                    AT_rr(:,m+1)=Lm*deltaCm;
                                    BT_rr(:,m+1)=Lm*Sm;                             
                                elseif volbaALFs==2
                                    if m==nmax
                                        T_rr=zeros(length_fiG,length_lambda);
                                    end

                                    T_rr=bsxfun(@times,T_rr,u)+(Lm*deltaCm*cos(m*lambda')+Lm*Sm*sin(m*lambda'));
                                end
    
                            elseif volbapar(i)==25 %Second radial derivative of disturbing potential

                                if m==nmax                               
                                    ampl_Wrr=((0:nmax)+1).*((0:nmax)+2);
                                end
                                
                                Lm=bsxfun(@times,Pnm(:,m_min:(nmax-m+1)),ampl_Wrr(:,(m_max+1):end));

                                if volbaALFs==1 || volbaALFs==3
                                    if m==nmax
                                        AWrr=zeros(length_fiG,nmax+1);
                                        BWrr=zeros(length_fiG,nmax+1);
                                    end

                                    AWrr(:,m+1)=Lm*Cm;
                                    BWrr(:,m+1)=Lm*Sm;                           
                                elseif volbaALFs==2
                                    if m==nmax
                                        Wrr=zeros(length_fiG,length_lambda);
                                    end

                                    Wrr=bsxfun(@times,Wrr,u)+(Lm*Cm*cos(m*lambda')+Lm*Sm*sin(m*lambda'));
                                end
                            end
                        end                                       
                    end

                    %Update of the progress bar
                    if GUI==1 %If working with the GUI
                        set(progressbar,'string',...
                            'Progress: Matrix multiplications...',...
                            'fontsize',8); drawnow;
                    end
                    
                    if GUI==0 && Status_bar==1 %If working without the GUI
                        fprintf('\n');
                        fprintf('Matrix multiplications...\n');
                    end

                    clear Lm dLm Pnm dPnm ddPnm Cm Sm C CEl CElm deltaC ...
                        deltaCm S u t q q2 index tu qu enm singdPnm singddPnm
                    
                    clear ampl_Trr ampl_Tfl ampl_Tzz amplH ampl_Vrr ampl_Vfl ...
                        ampl_Wr ampl_g_sa ampl_anomalia_sa ampl_Wrpor ...
                        ampl_porucha_sa ampl_T_rr ampl_Wrr
                    
                    if volbaALFs==1 || volbaALFs==3
                        cosla=cos((0:nmax)'*lambda');
                        sinla=sin((0:nmax)'*lambda');
                    end

                    %Computation of the normal gravity for eta, xi, Theta,
                    %Geoid undulation, zeta el, zeta
                    if any(volbapar==2) || any(volbapar==3) || any(volbapar==4) || any(volbapar==10) || any(volbapar==22) || any(volbapar==23)                        
                        bEl=aEl*sqrt(1-eEl^2);
                        EEl=sqrt(aEl^2-bEl^2);

                        %Computation of ellipsoidal harmonic coordinates
                        [X,Y,Z]=ell2cart(fi,0*zeros(length_fi,1),h,[aEl eEl]);

                        ugama=sqrt((X.*X+Y.*Y+Z.*Z-EEl.^2).*(1./2+1./2.*sqrt(1+((4.*EEl.^2.*Z.*Z)./(X.*X+Y.*Y+Z.*Z-EEl.^2).^2))));                        
                        betagama=atan((Z.*sqrt(ugama.^2+EEl.^2))./(ugama.*sqrt(X.*X+Y.*Y)));
                        clear X Y Z
                        
                        wgama=sqrt((ugama.^2+EEl^2*sin(betagama).^2)./(ugama.^2+EEl^2));
                        qgama=1/2*((1+(3*ugama.^2)./EEl^2).*atan(EEl./ugama)-3*ugama./EEl);
                        qgama_=3*(1+(ugama.^2)./EEl^2).*(1-ugama./EEl.*atan(EEl./ugama))-1;
                        qgama0=1/2*((1+(3*bEl^2)/EEl^2)*atan(EEl/bEl)-3*bEl/EEl);
                        
                        gamau=-1./wgama.*(GMEl./(ugama.^2+EEl.^2)+(omegaEl.^2.*aEl.^2.*EEl)./(ugama.^2+EEl.^2).*(qgama_./qgama0).*(1./2.*sin(betagama).^2-1./6)-omegaEl.^2.*ugama.*cos(betagama).^2);
                        gamabeta=-1./wgama.*(-(omegaEl.^2.*aEl.^2.*qgama)./(sqrt(ugama.^2+EEl.^2).*qgama0)+omegaEl.^2.*sqrt(ugama.^2+EEl.^2)).*sin(betagama).*cos(betagama);
                        
                        clear ugama betagama wgama qgama qgama_ qgama0
                        
                        gamaP=sqrt(gamau.^2+gamabeta.^2);
                        
                        clear gamau gamabeta
                    end

                    %% Final computation of functionals of the geopotential             
                    for i=1:pocetpar
                        if volbapar(i)==1                
                        elseif volbapar(i)==2 %Deflection of the vertical eta
                            if volbaALFs==1 || volbaALFs==3
                                eta=-Aeta*sinla+Beta*cosla;                                
                                clear Aeta Beta
                            elseif volbaALFs==2
                                eta=eta*1e280;
                            end

                            Pg=bsxfun(@times,-GM./(r.^2.*gamaP.*cos(fiG)),(eta))*(180/pi)*3600;
                            Pg(fi==pi/2 | fi==-pi/2,:)=0;
                            Pg=Pg(:);

                            clear eta
                        elseif volbapar(i)==3 %Deflection of the vertical xi
                            if volbaALFs==1 || volbaALFs==3
                                ksi=Aksi*cosla+Bksi*sinla;
                                clear Aksi Bksi
                            elseif volbaALFs==2
                                ksi=ksi*1e280;
                            end

                            Pg=bsxfun(@times,-GM./(r.^2.*gamaP),(ksi))*(180/pi)*3600;
                            Pg=Pg(:);

                            clear ksi
                        elseif volbapar(i)==4 %Deflection of the vertical Theta
                            if volbaALFs==1 || volbaALFs==3
                                Teta=-ATeta*sinla+BTeta*cosla;
                                clear ATeta BTeta
                                Tksi=ATksi*cosla+BTksi*sinla;
                                clear ATksi BTksi
                            elseif volbaALFs==2
                                Teta=Teta*1e280;
                                Tksi=Tksi*1e280;
                            end

                            Teta=bsxfun(@times,-GM./(r.^2.*gamaP.*cos(fiG)),(Teta))*(180/pi)*3600;
                            Teta(fi==pi/2 | fi==-pi/2,:)=0;
                            Tksi=bsxfun(@times,-GM./(r.^2.*gamaP),(Tksi))*(180/pi)*3600;
                            Teta=Teta(:);
                            Tksi=Tksi(:);
                            Talfa=atan2(Teta,Tksi);
                            Talfa(Talfa<0)=Talfa(Talfa<0)+2*pi;
   
                            Pg=[sqrt(Teta.^2+Tksi.^2) 180/pi*(Talfa)];

                            clear Teta Tksi Talfa
                        elseif volbapar(i)==5 %Disturbing potential                        
                            if volbaALFs==1 || volbaALFs==3
                                T=AT*cosla+BT*sinla;
                                clear AT BT
                            elseif volbaALFs==2
                                T=T*1e280;
                            end

                            Pg=bsxfun(@times,GM./r,T);
                            Pg=Pg(:);

                            clear T
                        elseif volbapar(i)==6 %Disturbing tensor Trr_Tpp_Tll
                            
                            clear Lmff Lmll                            
                            
                            if volbaALFs==1 || volbaALFs==3
                                Trr=ATrr*cosla+BTrr*sinla;
                                clear ATrr BTrr
                                Tff=ATff*cosla+BTff*sinla;
                                clear ATff BTff
                                
                                ATll(fi>pi/180*(89.9) | fi<pi/180*(-89.9),:)=0;
                                BTll(fi>pi/180*(89.9) | fi<pi/180*(-89.9),:)=0;
                                Tll=ATll*cosla+BTll*sinla;
                                clear ATll BTll
                            elseif volbaALFs==2
                                Trr=Trr*1e280;
                                Tff=Tff*1e280;
                                
                                Tll(fi>pi/180*(89.9) | fi<pi/180*(-89.9),:)=0; %#ok<*AGROW>
                                Tll=Tll*1e280;
                            end

                            Trr=bsxfun(@times,GM./r.^3,Trr)*10^9;
                            Tff=bsxfun(@times,GM./r.^3,Tff)*10^9;
                            Tll=bsxfun(@times,GM./(r.^3.*cos(fiG).^2),Tll)*10^9;
                            
                            Pg=Trr(:);
                            clear Trr
                            Pg=[Pg Tff(:)];
                            clear Tff
                            Pg=[Pg -Tll(:)];
                            clear Tll  
                            
                        elseif volbapar(i)==7 %Disturbing tensor Trp_Trl_Tpl
                            
                            clear Lmrf Lmrl Lmfl
                            
                            if volbaALFs==1 || volbaALFs==3
                                Trf=ATrf*cosla+BTrf*sinla;
                                clear ATrf BTrf
                                Trl=-ATrl*sinla+BTrl*cosla;
                                clear ATrl BTrl
                                
                                ATfl(fi>pi/180*(89.5) | fi<pi/180*(-89.5),:)=0;
                                BTfl(fi>pi/180*(89.5) | fi<pi/180*(-89.5),:)=0;
                                Tfl=-ATfl*sinla+BTfl*cosla;
                                clear ATfl BTfl
                            elseif volbaALFs==2
                                Trf=Trf*1e280;
                                Trl=Trl*1e280;
                                
                                Tfl(fi>pi/180*(89.5) | fi<pi/180*(-89.5),:)=0;
                                Tfl=Tfl*1e280;
                            end

                            Trf=bsxfun(@times,GM./r.^3,Trf)*10^9;
                            Trl=bsxfun(@times,GM./(r.^3.*cos(fiG)),Trl)*10^9;
                            Tfl=bsxfun(@times,GM./(r.^3.*cos(fiG)),Tfl)*10^9;
                            
                            Pg=-Trf(:);
                            clear Trf
                            Pg=[Pg -Trl(:)];
                            clear Trl
                            Pg=[Pg Tfl(:)];
                            clear Tfl  
                            
                        elseif volbapar(i)==8 %Disturbing tensor Txx_Tyy_Tzz
                            
                            clear anm bnm cnm LmTxx1 LmTxx2 LmTxx3 ...
                                LmTyy1 ampl_Tzz
                            
                            if volbaALFs==1 || volbaALFs==3
                                Tzz=ATzz*cosla+BTzz*sinla;
                                clear ATzz BTzz                                
                                
                                %Txx
                                Txx1=ATxx1*cosla+BTxx1*sinla;
                                clear ATxx1 BTxx1
                                
                                Txx2=ATxx2*[cosla(3:end,:);zeros(2,length(cosla(1,:)))]+BTxx2*[sinla(3:end,:);zeros(2,length(cosla(1,:)))];
                                clear ATxx2 BTxx2
                                
                                Txx=Txx1+Txx2;

                                %Tyy
                                Tyy1=ATyy1*cosla+BTyy1*sinla;
                                clear ATyy1 BTyy1
                                
                                Tyy=Tyy1+Txx2;
                                clear Txx1 Txx2 Tyy1
                                                         
                                Txx3=ATxx3*[zeros(2,length(cosla(1,:)));cosla(1:(end-2),:)]+BTxx3*[zeros(2,length(sinla(1,:)));sinla(1:(end-2),:)];
                                clear ATxx3 BTxx3
                                
                                Txx=Txx+Txx3;                                                                                                                       
                                Tyy=Tyy+Txx3;
                                clear Txx3
                                
                            elseif volbaALFs==2
                            end
                                                    
                            Txx=bsxfun(@times,GM./r.^3,Txx)*10^9;                           
                            Tyy=bsxfun(@times,GM./r.^3,Tyy)*10^9;
                            Tzz=bsxfun(@times,GM./r.^3,Tzz)*10^9; 
                                    
                            Pg=Txx(:);
                            clear Txx
                            Pg=[Pg -Tyy(:)];
                            clear Tyy
                            Pg=[Pg Tzz(:)];
                            clear Tzz

                        elseif volbapar(i)==9 %Disturbing tensor Txy_Txz_Tyz
                            
                            clear gamanm betanm gnm hnm dnm minm ninm ...
                                LmTxz1 LmTxz2 LmTxy1 LmTxy2 LmTxy3 ...
                                LmTyz1 LmTyz2
                            
                            if volbaALFs==1 || volbaALFs==3

                                %Txz
                                Txz1=ATxz1*[cosla(2:end,:);zeros(1,length(cosla(1,:)))]+BTxz1*[sinla(2:end,:);zeros(1,length(cosla(1,:)))];
                                clear ATxz1 BTxz1
                                
                                Txz2=ATxz2*[zeros(1,length(cosla(1,:)));cosla(1:(end-1),:)]+BTxz2*[zeros(1,length(sinla(1,:)));sinla(1:(end-1),:)];
                                clear ATxz2 BTxz2
                                
                                Txz=Txz1+Txz2;
                                clear Txz1 Txz2
                                
                                %Txy
                                Txy1=-ATxy1*[sinla(3:end,:);zeros(2,length(sinla(1,:)))]+BTxy1*[cosla(3:end,:);zeros(2,length(cosla(1,:)))];
                                clear ATxy1 BTxy1
                                
                                Txy2=-ATxy2*sinla+BTxy2*cosla;
                                clear ATxy2 BTxy2
                                
                                Txy3=-ATxy3*[zeros(2,length(sinla(1,:)));sinla(1:(end-2),:)]+BTxy3*[zeros(2,length(cosla(1,:)));cosla(1:(end-2),:)];
                                clear ATxy3 BTxy3

                                Txy=Txy1+Txy2+Txy3;
                                clear Txy1 Txy2 Txy3
                                
                                %Tyz
                                Tyz1=-ATyz1*[sinla(2:end,:);zeros(1,length(sinla(1,:)))]+BTyz1*[cosla(2:end,:);zeros(1,length(cosla(1,:)))];
                                clear ATyz1 BTyz1
           
                                Tyz2=-ATyz2*[zeros(1,length(sinla(1,:)));sinla(1:(end-1),:)]+BTyz2*[zeros(1,length(cosla(1,:)));cosla(1:(end-1),:)];
                                clear ATyz2 BTyz2
 
                                Tyz=Tyz1+Tyz2;
                                clear Tyz1 Tyz2
                            elseif volbaALFs==2
                            end

                            Txz=bsxfun(@times,GM./r.^3,Txz)*10^9;
                            Txy=bsxfun(@times,GM./r.^3,Txy)*10^9;
                            Tyz=bsxfun(@times,GM./r.^3,Tyz)*10^9;

                            Pg=Txy(:);
                            clear Txy
                            Pg=[Pg Txz(:)];
                            clear Txz
                            Pg=[Pg Tyz(:)];
                            clear Tyz

                        elseif volbapar(i)==10 %Geoid undulation
                            
                            clear HC HCm HS HSm LmH
                            
                            if volbaALFs==1 || volbaALFs==3
                                H=AH*cosla+BH*sinla;
                                clear AH BH
                                N1c=AN1c*cosla+BN1c*sinla;
                                clear AN1c BN1c
                            elseif volbaALFs==2
                                H=H*1e280;
                                N1c=N1c*1e280;
                            end                            
                            
                            N1c=bsxfun(@times,GM./(r.*gamaP),N1c);
                            H(H<0)=H(H<0)*0; %H is set to zero in the areas of oceans and seas
                           
                            G=6.67259*10^-11; %Newtonian gravitational constant (Moritz, 2000, Geodetic reference system 1980)
                            ro=2670; %Density of the crust

                            Pg=N1c-bsxfun(@times,(2*pi*G*ro*H.^2),1./gamaP);
                            Pg=Pg(:);
                            
                            clear H N1c
                        elseif volbapar(i)==11 %Gravitational potential
                            if volbaALFs==1 || volbaALFs==3
                                V=AV*cosla+BV*sinla;
                                clear AV BV
                            elseif volbaALFs==2
                                V=V*1e280;
                            end

                            Pg=bsxfun(@times,GM./r,V);
                            
                            Pg=Pg(:);

                            clear V 
                        elseif volbapar(i)==12 %Gravitational tensor Vrr_Vpp_Vll
                            
                            clear Lmff Lmll                            
                            
                            if volbaALFs==1 || volbaALFs==3
                                Vrr=AVrr*cosla+BVrr*sinla;
                                clear AVrr BVrr
                                Vff=AVff*cosla+BVff*sinla;
                                clear AVff BVff
                                
                                AVll(fi>pi/180*(89.9) | fi<pi/180*(-89.9),:)=0;
                                BVll(fi>pi/180*(89.9) | fi<pi/180*(-89.9),:)=0;
                                Vll=AVll*cosla+BVll*sinla;
                                clear AVll BVll
                            elseif volbaALFs==2
                                Vrr=Vrr*1e280;
                                Vff=Vff*1e280;
                                
                                Vll(fi>pi/180*(89.9) | fi<pi/180*(-89.9),:)=0;
                                Vll=Vll*1e280;
                            end

                            Vrr=bsxfun(@times,GM./r.^3,Vrr);
                            Vrr=Vrr*10^9;
                            
                            Vff=bsxfun(@times,GM./r.^3,Vff)*10^9;
                            Vll=bsxfun(@times,GM./(r.^3.*cos(fiG).^2),Vll)*10^9;
                            
                            Pg=Vrr(:);
                            clear Vrr
                            Pg=[Pg Vff(:)];
                            clear Vff
                            Pg=[Pg -Vll(:)];
                            clear Vll  
                            
                        elseif volbapar(i)==13 %Gravitational tensor Vrp_Vrl_Vpl
                            
                            clear Lmrf Lmrl Lmfl
                            
                            if volbaALFs==1 || volbaALFs==3
                                Vrf=AVrf*cosla+BVrf*sinla;
                                clear AVrf BVrf
                                Vrl=-AVrl*sinla+BVrl*cosla;
                                clear AVrl BVrl
                                
                                AVfl(fi>pi/180*(89.5) | fi<pi/180*(-89.5),:)=0;
                                BVfl(fi>pi/180*(89.5) | fi<pi/180*(-89.5),:)=0;
                                Vfl=-AVfl*sinla+BVfl*cosla;
                                clear AVfl BVfl
                            elseif volbaALFs==2
                                Vrf=Vrf*1e280;
                                Vrl=Vrl*1e280;
                                
                                Vfl(fi>pi/180*(89.5) | fi<pi/180*(-89.5),:)=0;
                                Vfl=Vfl*1e280;
                            end

                            Vrf=bsxfun(@times,GM./r.^3,Vrf)*10^9;
                            Vrl=bsxfun(@times,GM./(r.^3.*cos(fiG)),Vrl)*10^9;
                            Vfl=bsxfun(@times,GM./(r.^3.*cos(fiG)),Vfl)*10^9;
                            
                            Pg=-Vrf(:);
                            clear Vrf
                            Pg=[Pg -Vrl(:)];
                            clear Vrl
                            Pg=[Pg Vfl(:)];
                            clear Vfl  
                            
                        elseif volbapar(i)==14 %Gravitational tensor Vxx_Vyy_Vzz
                            
                            clear anm bnm cnm LmVxx1 LmVxx2 LmVxx3 ...
                                LmVyy1 ampl_Vzz
                            
                            if volbaALFs==1 || volbaALFs==3
                                Vzz=AVzz*cosla+BVzz*sinla;
                                clear AVzz BVzz                                
                                
                                %Vxx
                                Vxx1=AVxx1*cosla+BVxx1*sinla;
                                clear AVxx1 BVxx1
                                
                                Vxx2=AVxx2*[cosla(3:end,:);zeros(2,length(cosla(1,:)))]+BVxx2*[sinla(3:end,:);zeros(2,length(cosla(1,:)))];
                                clear AVxx2 BVxx2
                                
                                Vxx=Vxx1+Vxx2;
                                
                                %Vyy
                                Vyy1=AVyy1*cosla+BVyy1*sinla;
                                clear AVyy1 BVyy1
                                
                                Vyy=Vyy1+Vxx2;
                                clear Vxx1 Vxx2 Vyy1
                                                         
                                Vxx3=AVxx3*[zeros(2,length(cosla(1,:)));cosla(1:(end-2),:)]+BVxx3*[zeros(2,length(sinla(1,:)));sinla(1:(end-2),:)];
                                clear AVxx3 BVxx3
                                
                                Vxx=Vxx+Vxx3;                                                                                                                       
                                Vyy=Vyy+Vxx3;
                                clear Vxx3
                                
                            elseif volbaALFs==2
                            end

                            Vxx=bsxfun(@times,GM./r.^3,Vxx)*10^9;
                            Vyy=bsxfun(@times,GM./r.^3,Vyy)*10^9;
                            Vzz=bsxfun(@times,GM./r.^3,Vzz)*10^9;
                            
                            Pg=Vxx(:);
                            clear Vxx
                            Pg=[Pg -Vyy(:)];
                            clear Vyy
                            Pg=[Pg Vzz(:)];
                            clear Vzz
                            
                        elseif volbapar(i)==15 % Gravitational tensor Vxy_Vxz_Vyz
                            
                            clear gamanm betanm gnm hnm dnm minm ninm ...
                                LmVxz1 LmVxz2 LmVxy1 LmVxy2 LmVxy3 ...
                                LmVyz1 LmVyz2
                            
                            if volbaALFs==1 || volbaALFs==3

                                %Vxz
                                Vxz1=AVxz1*[cosla(2:end,:);zeros(1,length(cosla(1,:)))]+BVxz1*[sinla(2:end,:);zeros(1,length(cosla(1,:)))];
                                clear AVxz1 BVxz1
                                
                                Vxz2=AVxz2*[zeros(1,length(cosla(1,:)));cosla(1:(end-1),:)]+BVxz2*[zeros(1,length(sinla(1,:)));sinla(1:(end-1),:)];
                                clear AVxz2 BVxz2
                                
                                Vxz=Vxz1+Vxz2;
                                clear Vxz1 Vxz2
                                
                                %Vxy
                                Vxy1=-AVxy1*[sinla(3:end,:);zeros(2,length(sinla(1,:)))]+BVxy1*[cosla(3:end,:);zeros(2,length(cosla(1,:)))];
                                clear AVxy1 BVxy1
                                
                                Vxy2=-AVxy2*sinla+BVxy2*cosla;
                                clear AVxy2 BVxy2
                                
                                Vxy3=-AVxy3*[zeros(2,length(sinla(1,:)));sinla(1:(end-2),:)]+BVxy3*[zeros(2,length(cosla(1,:)));cosla(1:(end-2),:)];
                                clear AVxy3 BVxy3

                                Vxy=Vxy1+Vxy2+Vxy3;
                                clear Vxy1 Vxy2 Vxy3
                                
                                %Vyz
                                Vyz1=-AVyz1*[sinla(2:end,:);zeros(1,length(sinla(1,:)))]+BVyz1*[cosla(2:end,:);zeros(1,length(cosla(1,:)))];
                                clear AVyz1 BVyz1
                                
                                Vyz2=-AVyz2*[zeros(1,length(sinla(1,:)));sinla(1:(end-1),:)]+BVyz2*[zeros(1,length(cosla(1,:)));cosla(1:(end-1),:)];
                                clear AVyz2 BVyz2
 
                                Vyz=Vyz1+Vyz2;
                                clear Vyz1 Vyz2
                            elseif volbaALFs==2
                            end

                            Vxz=bsxfun(@times,GM./r.^3,Vxz)*10^9;
                            Vxy=bsxfun(@times,GM./r.^3,Vxy)*10^9;
                            Vyz=bsxfun(@times,GM./r.^3,Vyz)*10^9;

                            Pg=Vxy(:);
                            clear Vxy
                            Pg=[Pg Vxz(:)];
                            clear Vxz
                            Pg=[Pg Vyz(:)];
                            clear Vyz

                        elseif volbapar(i)==16 %Gravity vector gX_gY_gZ
                            clear LmWr LmWlambda LmWfi
                            
                            if volbaALFs==1 || volbaALFs==3
                                Wr=AWr*cosla+BWr*sinla;
                                clear AWr BWr
                                Wlambda=-AWlambda*sinla+BWlambda*cosla;
                                clear AWlambda BWlambda 
                                Wfi=AWfi*cosla+BWfi*sinla;
                                clear AWfi BWfi 
                            elseif volbaALFs==2
                                Wr=Wr*1e280;
                                Wlambda=Wlambda*1e280;
                                Wfi=Wfi*1e280;
                            end

                            Wr=bsxfun(@times,GM./r.^2,Wr);
                            Wr=bsxfun(@plus,-Wr,omegaEl^2.*r.*(cos(fiG).^2));
                            Wr=Wr(:);

                            Wlambda=bsxfun(@times,GM./r,Wlambda);
                            Wlambda=bsxfun(@times,-Wlambda,1./(r.*cos(fiG)));
                            Wlambda=Wlambda(:);

                            Wfi=bsxfun(@times,GM./r,Wfi);
                            Wfi=bsxfun(@plus,Wfi,-omegaEl^2.*(r.^2).*cos(fiG).*sin(fiG));
                            Wfi=bsxfun(@times,Wfi,1./r);
                            Wfi=Wfi(:);
                            
                            Pg=Wfi*10^5;
                            clear Wfi
                            Pg=[Pg Wlambda*10^5];
                            clear Wlambda
                            Pg=[Pg Wr*10^5];
                            clear Wr
                        elseif volbapar(i)==17 %Gravity sa
                            if volbaALFs==1 || volbaALFs==3
                                g_sa=Ag_sa*cosla+Bg_sa*sinla;
                                clear Ag_sa Bg_sa
                            elseif volbaALFs==2
                                g_sa=g_sa*1e280;
                            end

                            Pg=bsxfun(@times,GM./r.^2,g_sa);

                            Pg=sqrt(bsxfun(@plus,-Pg,omegaEl^2.*r.*(cos(fiG).^2)).^2)*10^5;                    
                            Pg=Pg(:);

                            clear g_sa
                        elseif volbapar(i)==18 %Gravity potential
                            if volbaALFs==1 || volbaALFs==3
                                W=AW*cosla+BW*sinla;
                                clear AW BW
                            elseif volbaALFs==2
                                W=W*1e280;
                            end

                            Pg=bsxfun(@times,GM./r,W);
                            clear W

                            Pg=bsxfun(@plus,Pg,1/2*omegaEl.^2.*r.^2.*cos(fiG).^2);
                            Pg=Pg(:);
                        elseif volbapar(i)==19 %Gravity anomaly sa
                            if volbaALFs==1 || volbaALFs==3
                                anomalia_sa=Aanomalia_sa*cosla+Banomalia_sa*sinla;
                                clear Aanomalia_sa Banomalia_sa
                            elseif volbaALFs==2
                                anomalia_sa=anomalia_sa*1e280;
                            end

                            Pg=bsxfun(@times,GM./r.^2,anomalia_sa)*10^5;
                            Pg=Pg(:);

                            clear anomalia_sa
                        elseif volbapar(i)==20 %Gravity disturbance
                            clear LmWrpor LmWlambdapor LmWfipor
                            
                            if volbaALFs==1 || volbaALFs==3
                                Wrpor=AWrpor*cosla+BWrpor*sinla;
                                clear AWrpor BWrpor
                                Wlambdapor=-AWlambdapor*sinla+BWlambdapor*cosla;
                                clear AWlambdapor BWlambdapor
                                Wfipor=AWfipor*cosla+BWfipor*sinla;
                                clear AWfipor BWfipor
                                Ur=AUr*cos(0*lambda');
                                clear AUr
                                Ufi=AUfi*cos(0*lambda');
                                clear AUfi
                            elseif volbaALFs==2
                                Wrpor=Wrpor*1e280;
                                Wlambdapor=Wlambdapor*1e280;
                                Wfipor=Wfipor*1e280;
                                Ur=Ur*1e280;
                                Ufi=Ufi*1e280;
                            end

                            Wrpor=bsxfun(@times,GM./r.^2,Wrpor);
                            Wrpor=sqrt(bsxfun(@plus,-Wrpor,omegaEl^2.*r.*(cos(fiG).^2)).^2);                    

                            Wlambdapor=bsxfun(@times,GM./r,Wlambdapor);
                            Wlambdapor=bsxfun(@times,Wlambdapor,1./(r.*cos(fiG)));

                            Wfipor=bsxfun(@times,GM./r,Wfipor);
                            Wfipor=bsxfun(@plus,Wfipor,-omegaEl^2.*(r.^2).*cos(fiG).*sin(fiG));
                            Wfipor=bsxfun(@times,Wfipor,1./r);

                            Ur=bsxfun(@times,GM./r.^2,Ur);
                            Ur=sqrt(bsxfun(@plus,-Ur,omegaEl^2.*r.*(cos(fiG).^2)).^2);                    

                            Ufi=bsxfun(@times,GM./r,Ufi);
                            Ufi=bsxfun(@plus,Ufi,-omegaEl^2.*(r.^2).*cos(fiG).*sin(fiG));
                            Ufi=bsxfun(@times,Ufi,1./r);

                            Pg=(sqrt(Wrpor.^2+Wlambdapor.^2+Wfipor.^2)-sqrt(Ur.^2+Ufi.^2))*10^5;
                            Pg=Pg(:);
                            
                            clear Wrpor Wlambdapor Wfipor Ur Ufi
                        elseif volbapar(i)==21 %Gravity disturbance sa
                            if volbaALFs==1 || volbaALFs==3
                                porucha_sa=Aporucha_sa*cosla+Bporucha_sa*sinla;
                                clear Aporucha_sa Bporucha_sa
                            elseif volbaALFs==2
                                porucha_sa=porucha_sa*1e280;
                            end

                            Pg=bsxfun(@times,GM./r.^2,porucha_sa)*10^5;
                            Pg=Pg(:);

                            clear porucha_sa 
                        elseif volbapar(i)==22 %Height anomaly ell
                            if volbaALFs==1 || volbaALFs==3
                                zetaEl=AzetaEl*cosla+BzetaEl*sinla;
                                clear AzetaEl BzetaEl
                            elseif volbaALFs==2
                                zetaEl=zetaEl*1e280;
                            end

                            Pg=bsxfun(@times,GM./(r.*gamaP),zetaEl);
                            Pg=Pg(:);

                            clear zetaEl 
                        elseif volbapar(i)==23 %Height anomaly
                            
                            clear HC HCm HS HSm LmH Lmdg
                            
                            if volbaALFs==1 || volbaALFs==3
                                zeta_H=AH_zeta*cosla+BH_zeta*sinla;
                                clear AH_zeta BH_zeta
                                zeta_N1c=AN1c_zeta*cosla+BN1c_zeta*sinla;
                                clear AN1c_zeta BN1c_zeta
                                zeta_dg=Azetadg*cosla+Bzetadg*sinla;
                                clear Azetadg Bzetadg
                                zeta_zetaEl=Azeta*cosla+Bzeta*sinla;
                                clear Azeta Bzeta  
                            elseif volbaALFs==2
                                zeta_H=zeta_H*1e280;
                                zeta_N1c=zeta_N1c*1e280;
                                zeta_dg=zeta_dg*1e280;
                                zeta_zetaEl=zeta_zetaEl*1e280;
                            end                                 
                            
                            zeta_N1c=bsxfun(@times,GM./(r.*gamaP),zeta_N1c);
                            zeta_H(zeta_H<0)=zeta_H(zeta_H<0)*0; %H is set to zero in the areas of oceans and seas
   
                            zeta_zetaEl=bsxfun(@times,GM./(r.*gamaP),zeta_zetaEl);
                            
                            zeta_dg=bsxfun(@times,GM./r.^2,zeta_dg);
                            
                            G=6.67259*10^-11; %Newtonian gravitational constant (Moritz, 2000, Geodetic reference system 1980)
                            ro=2670; %Density of the crust
                            
                            zeta_N=(zeta_N1c-bsxfun(@times,(2*pi*G*ro*zeta_H.^2),1./gamaP));

                            Pg=zeta_zetaEl-bsxfun(@times,zeta_dg.*(zeta_H+zeta_N),1./gamaP);
                            Pg=Pg(:);
                            
                            clear zeta_N1c zeta_H zeta_zetaEl zeta_dg zeta_N
                        elseif volbapar(i)==24 %Second radial derivative of disturbing potential
                            if volbaALFs==1 || volbaALFs==3
                                T_rr=AT_rr*cosla+BT_rr*sinla;
                                clear AT_rr BT_rr
                            elseif volbaALFs==2
                                T_rr=T_rr*1e280;
                            end

                            Pg=bsxfun(@times,GM./r.^3,T_rr)*10^9;                            
                            Pg=Pg(:);

                            clear T_rr
                        elseif volbapar(i)==25 %Second radial derivative of gravity potential
                            if volbaALFs==1 || volbaALFs==3
                                Wrr=AWrr*cosla+BWrr*sinla;
                                clear AWrr BWrr
                            elseif volbaALFs==2
                                Wrr=Wrr*1e280;
                            end

                            Pg=bsxfun(@times,GM./r.^3,Wrr);
                            clear Wrr
                            Pg=bsxfun(@plus,Pg,omegaEl^2.*cos(fiG).^2)*10^9;
                            Pg=Pg(:);
                        end

                        if i==1
                            P=Pg;
                            clear Pg
                        else
                            P=[P Pg];
                            if i==pocetpar
                                clear Pg
                            end
                        end
                    end

                    %Update of the progress bar
                    set(progressbar,'string','','fontsize',8); drawnow;
                    
                    clear sinla cosla r gamaP                                
                end                

                %% Computation of functionals in the point-wise mode
                %==============================================================

                %Identification of point-wise approach                         
                if volbadiskcheck==1                           
                    %Ellipsoidal coordinates
                    if GUI==1 %If working with the GUI
                        fi=str2num(get(findobj('tag','fi'),'string'))';
                        lambda=str2num(get(findobj('tag','lambda'),'string'))';
                        h=str2num(get(findobj('tag','hdisk'),'string'))';
                    elseif GUI==0 %If working without the GUI
                        fi=lat(:);
                        lambda=lon(:);
                        h=h2(:);
                        clear lat lon h2
                    end
                end
                
                %Identification of load data approach
                if volbaloadcheck==1

                    if isempty(loadname) %Error message, if the file with input points
                        %has not been imported
                        if GUI==1 %If working with the GUI
                            errordlg('Please input the data file containing coordinates of the computing points.',...
                                'Error in point type selection');
                        end
                        error('Please input the data file containing coordinates of the computing points.')
                    end
                    
                    if ~exist([loadadresar,loadname],'file') %Check whether the input data point file exists
                        if GUI==1 %Working with the GUI
                            errordlg('The entered data point file does not exist.',...
                                'Error in point type selection');
                        end
                        error('The entered data point file does not exist.')
                    end
                    if strcmp(loadname(end-3:end),'.mat') %Loading MAT file
                        Import=load([loadadresar,loadname]);
                        Import=struct2cell(Import);
                        Import=cell2mat(Import);
                        Import=Import(:,1:3);
                    else
                        Import=load([loadadresar,loadname]);
                    end

                    [rows_Import,cols_Import]=size(Import);
                    
                    if cols_Import<3
                        if GUI==1 %If working with the GUI
                            errordlg('The input data file containing coordinates of the computational points must have three columns.',...
                                'Error in point type selection');
                        end
                        error('The input data file containing coordinates of the computational points must have three columns.')
                    end
                    
                    clear rows_Import cols_Import
                    
                    %Ellipsoidal coordinates
                    fi=Import(:,1);
                    lambda=Import(:,2);
                    h=Import(:,3);

                    volbadiskcheck=1;
                end
                
                length_fi=length(fi);
                length_lambda=length(lambda);
                %=============================================================
                if volbadiskcheck==1                             

                    %Error message for input ellipsoidal coordinates
                    if isempty(fi)
                        if GUI==1
                            errordlg('The "Latitude" array cannot be empty.',...
                                'Error in point type selection');
                            error('The "Latitude" array cannot be empty.');
                        elseif GUI==0
                            error('The "lat" variable cannot be empty.');
                        end
                    end
                    if isempty(lambda)
                        if GUI==1
                            errordlg('The "Longitude" array cannot be empty.',...
                                'Error in point type selection');
                            error('The "Longitude" array cannot be empty.');
                        elseif GUI==0
                            error('The "lon" variable cannot be empty.');
                        end
                    end
                    if isempty(h)
                        if GUI==1
                            errordlg('The "Ellipsoidal height/Spherical radius" array cannot be empty.',...
                                'Error in point type selection');
                            error('The "Ellipsoidal height/Spherical radius" array cannot be empty.');
                        elseif GUI==0
                            error('The "h2" variable cannot be empty.');
                        end
                    end

                    if length_fi~=length_lambda || length_fi~=length(h) || length_lambda~=length(h)
                        if GUI==1 %If working with the GUI
                            errordlg('Coordinates dimensions are not consistent.',...
                                'Error in point type selection')
                        end
                        error('Coordinates dimensions are not consistent.')     
                    end               

                    if any(fi>90) || any(fi<-90)
                        if GUI==1 %If working with the GUI
                            errordlg('The values of latitude must be within the interval <-90 deg , 90 deg>.',...
                                'Error in point type selection');
                            error('The values of latitude must be within the interval <-90 deg, 90 deg>.');
                        elseif GUI==0
                            error('The values of the "lat" variable must be within the interval <-90 deg, 90 deg>.');
                        end
                    end
                    if any(lambda>360) || any(lambda<-180)
                        if GUI==1 %If working with the GUI
                            errordlg('The values of longitudes must be within the interval <-180 deg, 180 deg> or <0 deg, 360 deg>.',...
                                'Error in point type selection');
                            error('The values of longitudes must be within the interval <-180 deg, 180 deg> or <0 deg, 360 deg>.');
                        elseif GUI==0
                            error('The values of the "lon" variable must be within the interval <-180 deg, 180 deg> or <0 deg, 360 deg>.');
                        end
                    end

                    fi=pi/180*(fi(:));
                    lambda=pi/180*(lambda(:));               

                    if coord==1 %Entered spherical coordinates                       
                        %Spherical radius
                        r=h;

                        %Spherical latitude
                        fiG=fi;  
                        
                        %Transform spherical latitude into the ellipsoidal
                        %latitude
                        [X,Y,Z]=sph2cart(lambda,fiG,r);
                        [fi,lambda_del,h]=cart2ell(X,Y,Z,[aEl eEl]);

                        clear X Y Z lambda_del
                    elseif coord==0 %Entered ellipsoidal coordinates
                        %Trasformation of (fi, lambda, h) into (X, Y, Z)
                        [X,Y,Z]=ell2cart(fi,lambda,h,[aEl eEl]);  
                        r=sqrt(X.*X+Y.*Y+Z.*Z); %Radius

                        %Spherical latitude
                        fiG=atan(Z./sqrt(X.*X+Y.*Y)); 

                        clear X Y Z 
                    end

                    %Computation of the coefficients C0,0; C2,0; ...; C20,0
                    %of the selected ellipsoid
                    CEl=zeros(length(C),1);
                    for n=0:10
                        CEl(2*n==stupen & rad==0,1)=((-1)^n*(3*eEl^(2*n))/((2*n+1)*(2*n+3)*sqrt(4*n+1))*(1-n-5^(3/2)*n*CEl_20/eEl^2)).*(aEl./R).^(2*n).*(GMEl/GM);
                    end                                

                    if any(volbapar==11) || any(volbapar==12) || any(volbapar==13) || any(volbapar==14) || any(volbapar==15) || any(volbapar==16) || any(volbapar==17) || any(volbapar==18) || any(volbapar==20) || any(volbapar==25)
                        grav=1;
                    else
                        grav=0;
                    end

                    if  any(volbapar==2) || any(volbapar==3) || any(volbapar==4) || any(volbapar==5) || any(volbapar==6) || any(volbapar==7) || any(volbapar==8) || any(volbapar==9) || any(volbapar==10) || any(volbapar==19) || any(volbapar==21) || any(volbapar==22) || any(volbapar==23) || any(volbapar==24)
                        por=1;
                        deltaC=C-CEl;
                    else
                        por=0;
                    end

                    if any(volbapar==20)
                        normal=1;
                    else
                        normal=0;
                    end

                    clear GGM stupen rad
                    if normal==0
                        clear CEl
                    end                                      

                    %Initialization
                    eta=0; ksi=0; Teta=0; Tksi=0; T=0; T_rr=0; Trr=0; Trfi=0; Trl=0; Tfifi=0;
                    Tfil=0; Tll=0; N=0; V=0; Vrr=0; Vrfi=0; Vrl=0; Vfifi=0;
                    Vfil=0; Vll=0; g=0; g_sa=0; W=0; anomalia_sa=0; porucha=0; 
                    porucha_sa=0; zetaEl=0; zeta=0; Wrr=0; Wr=0; Wfi=0; Wlambda=0;
                    Ur=0; Ufi=0; N1c=0; N2c=0; H=0; zeta_N1c=0; zeta_H=0; zeta_dg=0;
                    etaH=0; ksiH=0; TetaH=0; TksiH=0; TH=0; T_rrH=0; TrrH=0; TrfiH=0; TrlH=0; TfifiH=0;
                    TfilH=0; TllH=0; NH=0; VH=0; VrrH=0; VrfiH=0; VrlH=0; VfifiH=0;
                    VfilH=0; VllH=0; gH=0; g_saH=0; WH=0; anomalia_saH=0; poruchaH=0; 
                    porucha_saH=0; zetaElH=0; zetaH=0; WrrH=0; WrH=0; WfiH=0; WlambdaH=0;
                    UrH=0; UfiH=0; N1cH=0; N2cH=0; HH=0; zeta_N1cH=0; zeta_HH=0; zeta_dgH=0;
                    Tzz=0; Txx=0; Tyy=0; Txy=0; Txz=0; Tyz=0;
                    TzzH=0; TxxH=0; TyyH=0; TxyH=0; TxzH=0; TyzH=0;
                    Vzz=0; Vxx=0; Vyy=0; Vxy=0; Vxz=0; Vyz=0;
                    VzzH=0; VxxH=0; VyyH=0; VxyH=0; VxzH=0; VyzH=0;

                    %Initialization of the matrices and vectors for the computation of fnALFs
                    length_fiG=length(fiG);
                    Pnm=zeros(length_fiG,nmax+1);
                    q=(R./r);
                    q2=(R./r).^2;
                    u=cos(fiG);
                    t=sin(fiG);
                    %Initialization for extended-range arithmetic approach
                    if volbaALFs==3
                        
                        bit=mexext; %Bit version of Matlab
                        bit=bit(end-1:end);
                        bit=str2double(bit);
                        if bit==32
                            bit=32;
                        elseif bit==64
                            bit=64;
                        else
                            bit=64;
                        end
                        
                        nmax23=nmax*2+3;
                        rr=zeros(nmax23,1); ri=rr;
                        dd=zeros(nmax,1); am=dd; bm=am;

                        m1=1:nmax23;
                        rr(m1)=sqrt(m1);
                        ri(m1)=1./rr;
                        m2=1:nmax;
                        dd(m2)=rr(2*m2+3).*ri(2*m2+2);

                        IND=960;
                        BIG=2^IND;
                        BIGI=2^(-IND);
                        BIGS=2^(IND/2);
                        BIGSI=2^(-IND/2);
                        ROOT3=1.732050807568877;
                        
                        if bit==32
                            pm=am;
                            ps1=zeros(length_fiG,nmax); 
                            ips1=ps1;
                            x=ROOT3*u.*q;
                            ix=zeros(size(x));
                            ps1(:,1)=x;
                            ips1(:,1)=ix;
                            for m3=2:nmax
                                x=(dd(m3-1)*u).*x.*q;
                                y=abs(x);
                                iy=y>=BIGS;
                                if any(iy)
                                    x(iy)=x(iy)*BIGI;
                                    ix(iy)=ix(iy)+1;
                                end
                                iy=y<BIGSI;
                                if any(iy)
                                    x(iy)=x(iy)*BIG;
                                    ix(iy)=ix(iy)-1;
                                end
                                ps1(:,m3)=x;
                                ips1(:,m3)=ix;
                            end
                        elseif bit==64
                            tq=t.*q;
                            temp1=zeros(length_fiG,1);
                            temp2=ones(length_fiG,1);
                            temp3=temp2;
                            temp4=temp1;
                            temp5=temp1+BIGI;
                            ps1b=zeros(length_fiG,nmax); 
                            ips1b=ps1b;
                            xb=ROOT3*u.*q;
                            ixb=zeros(size(xb));
                            ps1b(:,1)=xb;
                            ips1b(:,1)=ixb;
                            for m3=2:nmax
                                xb=(dd(m3-1)*u).*xb.*q;
                                yb=abs(xb);
                                iyb=yb>=BIGS;
                                if any(iyb)
                                    xb(iyb)=xb(iyb)*BIGI;
                                    ixb(iyb)=ixb(iyb)+1;
                                end
                                iyb=yb<BIGSI;
                                if any(iyb)
                                    xb(iyb)=xb(iyb)*BIG;
                                    ixb(iyb)=ixb(iyb)-1;
                                end
                                ps1b(:,m3)=xb;
                                ips1b(:,m3)=ixb;
                            end
                        end
                        
                        clear dd
                    end

                    %If geoid/quasigeoid is to be computed
                    if any(volbapar==10) || any(volbapar==23)
                        geoid=1;
                        if any(h~=0)
                            if GUI==1 %If working with the GUI
                                errordlg('To compute Geoid_undulation or Height_anomaly, the ellipsoidal height must be set to zero.',...
                                    'Error in point type selection');
                            end
                            error('To compute Geoid_undulation or Height_anomaly, the ellipsoidal height must be set to zero.');
                        end
                    else
                        geoid=0;
                    end

                    %Initialization of the matrices and vectors for the 
                    %computation of the first-order derivatives of fnALFs
                    if any(volbapar==3) || any(volbapar==4) || any(volbapar==6) || any(volbapar==7) || any(volbapar==12) || any(volbapar==13) || any(volbapar==16)  || any(volbapar==20)
                        dALFs=1;
                        dPnm=zeros(length_fi,nmax+1);
                        qu=q./u;
                        tu=t./u;
                        
                        %Treatment of the dPnm singularity
                        singdPnm=fi==pi/2 | fi==-pi/2;
                    else
                        dALFs=0;
                    end   
                    
                    %Initialization of the matrices and vectors for the 
                    %computation of the second-order derivatives of fnALFs
                    if any(volbapar==6) || any(volbapar==12)
                        ddALFs=1;
                        ddPnm=zeros(length_fi,nmax+1);
                        
                        %Treatment of the ddPnm singularity
                        singddPnm=fi==pi/2 | fi==-pi/2;
                    else
                        ddALFs=0;
                    end   

                    %Status line
                    progressbar=findobj('tag','hlasky');
                    
                    if GUI==0 && Status_bar==1 %If working without the GUI
                        fprintf('Progress: m = ')
                    end

                    
                    %% Summation over m
                    for m=nmax:-1:0

                        %Update of the progress bar
                        if GUI==1 && rem(m,10)==0 %If working with the GUI
                            set(progressbar,'string',...
                                sprintf('Progress: m = %5.0d',m),...
                                'fontsize',8); drawnow;
                        end
                        
                        if GUI==0 && Status_bar==1 %If working without the GUI
                            if rem(m,10)==0
                                fprintf('%d,',m)
                            end
                        end

                        m_min=1+max([0 nminGGM-m]);
                        
                        %Selection of the spherical harmonic coefficients of order m
                        %======================================================
                        if grav==1 %C's spherical harmonic coefficients for the functionals without the normal gravity field
                            Cm=C(index((m+m_min-nminGGM):end)+m);
                        end

                        if por==1 %C's spherical harmonic coefficients for the functionals with the disturbing field
                            deltaCm=deltaC(index((m+m_min-nminGGM):end)+m);
                        end
                                
                        if normal==1 %C's spherical harmonic coefficients for the functionals with the normal field
                            if m==0
                                CElm=CEl(index((m+m_min-nminGGM):end)+m);
                            end
                        end
                            
                        if geoid==1
                            HCm=HC(index((m+m_min-nminGGM):end)+m); 
                            HSm=HS(index((m+m_min-nminGGM):end)+m);
                        end

                        Sm=S(index((m+m_min-nminGGM):end)+m);
                        %======================================================

                        %% Computation of modified fnALFs
                        if volbaALFs==1 %Standard forward column method
                            if m==0
                                Pnm(:,1)=1;
                            elseif m==1                   
                                Pnm(:,1)=sqrt(3)*u.*q;  
                            elseif m>1                            
                                i=2*(2:m);
                                i1=sqrt((i+ones(size(i)))./i);
                                Pnm(:,1)=u.^m*sqrt(3)*prod(i1).*q.^m;
                            end

                            if m==nmax
                            elseif m<=(nmax-1)
                                n=m+1;
                                anm=sqrt((2*n-1)*(2*n+1)/((n-m)*(n+m)));
                                Pnm(:,2)=anm*t.*Pnm(:,1).*q;
                            end

                            if m<(nmax-1)
                                j=3;
                                for n=m+2:nmax                            
                                    anm=sqrt((2*n-1)*(2*n+1)/((n-m)*(n+m)));
                                    bnm=sqrt((2*n+1)*(n+m-1)*(n-m-1)/((n-m)*(n+m)*(2*n-3)));
                                    Pnm(:,j)=anm*t.*Pnm(:,j-1).*q-bnm*Pnm(:,j-2).*q2;
                                    j=j+1;
                                end
                            end

                        elseif volbaALFs==2 %Modified forward column method
                            if m==0
                                Pnm(:,1)=1e-280;
                            elseif m==1

                                Pnm(:,1)=sqrt(3)*q*1e-280;  
                            elseif m>1                            
                                i=2*(2:m);
                                i1=sqrt((i+ones(size(i)))./i);
                                Pnm(:,1)=sqrt(3)*prod(i1)*(q.^m)*1e-280;
                            end

                            if m==nmax
                            elseif m<=(nmax-1)
                                n=m+1;
                                anm=sqrt((2*n-1)*(2*n+1)/((n-m)*(n+m)));
                                Pnm(:,2)=anm*t.*Pnm(:,1).*q;
                            end

                            if m<(nmax-1)
                                j=3;
                                for n=m+2:nmax                            
                                    anm=sqrt((2*n-1)*(2*n+1)/((n-m)*(n+m)));
                                    bnm=sqrt((2*n+1)*(n+m-1)*(n-m-1)/((n-m)*(n+m)*(2*n-3)));
                                    Pnm(:,j)=anm*t.*Pnm(:,j-1).*q-bnm*Pnm(:,j-2).*q2;
                                    j=j+1;
                                end
                            end

                        elseif volbaALFs==3 %Extended-range arithmetic 
                            if bit==32 %32 bit version of Matlab
                                am(m+1)=rr(2*m+3);
                                for n=m+2:nmax
                                    w=rr(2*n+1)*ri(n-m)*ri(n+m);
                                    am(n)=rr(2*n-1)*w;
                                    bm(n)=rr(n-m-1)*rr(n+m-1)*ri(2*n-3)*w;
                                end

                                if m~=0
                                    for i=1:length_fiG                                   
                                        x=ps1(i,m);
                                        ix=ips1(i,m);
                                        if(ix==0)
                                            pm(m)=x;
                                        elseif (ix<-1)    
                                            pm(m)=0;         
                                        elseif (ix<0)
                                            pm(m)=x*BIGI;
                                        else
                                            pm(m)=x*BIG;
                                        end

                                        if(m>=nmax)
                                            Pnm(i,1:(nmax-m+1))=pm(m:end);
                                            continue;
                                        end
                                        y=x;
                                        iy=ix;
                                        x=(am(m+1)*t(i)*q(i))*y;
                                        ix=iy;
                                        w=abs(x);

                                        if(w>=BIGS)
                                            x=x*BIGI;
                                            ix=ix+1;
                                        elseif (w<BIGSI)
                                            x=x*BIG;
                                            ix=ix-1;
                                        end
                                        if(ix==0)
                                            pm(m+1)=x;
                                        elseif (ix<-1)    
                                            pm(m+1)=0.;        
                                        elseif (ix<0)
                                            pm(m+1)=x*BIGI;
                                        else
                                            pm(m+1)=x*BIG;
                                        end

                                        for n=m+2:nmax                                       
                                            id=ix-iy;
                                            if(id==0)
                                                zz=(am(n)*t(i)*q(i))*x-bm(n)*q2(i)*y;
                                                iz=ix;
                                            elseif (id==1)
                                                zz=(am(n)*t(i)*q(i))*x-bm(n)*q2(i)*(y*BIGI);
                                                iz=ix;
                                            elseif (id==-1)
                                                zz=(am(n)*t(i)*q(i))*(x*BIGI)-bm(n)*q2(i)*y;
                                                iz=iy;
                                            elseif (id>1)
                                                zz=(am(n)*t(i)*q(i))*x;
                                                iz=ix;
                                            else
                                                zz=-bm(n)*q2(i)*y;
                                                iz=iy;
                                            end

                                            w=abs(zz);
                                            if(w>=BIGS)
                                                zz=zz*BIGI;
                                                iz=iz+1;
                                            elseif (w<BIGSI)
                                                zz=zz*BIG;
                                                iz=iz-1;
                                            end

                                            if(iz==0)
                                                pm(n)=zz;
                                            elseif (iz<-1)   
                                                pm(n)=0.;          
                                            elseif (iz<0)
                                                pm(n)=zz*BIGI;
                                            else
                                                pm(n)=zz*BIG;
                                            end

                                            y=x;
                                            iy=ix;
                                            x=zz;
                                            ix=iz; 
                                        end

                                        Pnm(i,1:(nmax-m+1))=pm(m:end);         
                                    end

                                elseif m==0
                                    Pnm(:,1)=1;
                                    Pnm(:,2)=sqrt(3)*t.*q;

                                    for i=2:nmax
                                        Pnm(:,i+1)=Pnm(:,i).*sqrt((2*i+1)*(2*i-1))./i.*t.*q-q2.*Pnm(:,i-1).*(i-1).*sqrt(2.*i+1)./(i.*sqrt(2.*i-3));
                                    end

                                    clear rr ri dd am bm pm ps1 ips1 m1 m2 ...
                                        dd ix x y iy w iz z
                                end
                            else %64 bit version of Matlab
                                am(m+1)=rr(2*m+3);
                                for n=m+2:nmax
                                    w=rr(2*n+1)*ri(n-m)*ri(n+m);
                                    am(n)=rr(2*n-1)*w;
                                    bm(n)=rr(n-m-1)*rr(n+m-1)*ri(2*n-3)*w;
                                end
                                                              
                                if m==0 %Zonal modified fnALFs
                                    Pnm(:,1)=1;
                                    Pnm(:,2)=sqrt(3)*t.*q;

                                    for i=2:nmax
                                        Pnm(:,i+1)=Pnm(:,i).*sqrt((2*i+1)*(2*i-1))./i.*t.*q-q2.*Pnm(:,i-1).*(i-1).*sqrt(2.*i+1)./(i.*sqrt(2.*i-3));
                                    end

                                    clear rr ri am bm ps1 ips1 m1 m2 dd ...
                                        ix x y iy w iz zz pmx pm0 pmxBIGI ...
                                        pmxBIG w wBIGS wBIGSI pm1x pm10 ...
                                        pm1xBIGI pm1xBIG id id0 id1 id_1 ...
                                        idv1 idm1 iz0 izm_1 izm0 izv0

                                elseif m~=0 %Non-zonal modified fnALFs
                                    xb=ps1b(:,m);
                                    ixb=ips1b(:,m);
                                                                       
                                    temp5(ixb==0)=1;
                                    temp5(ixb<-1)=0;
                                    %temp5(izb>=-1 & izb<0)=BIGI;
                                    %The condition "izb>=-1 & izb<0"
                                    %is useless, as "izb" is already
                                    %initialized as "izb=BIGI".
                                    temp5(ixb>0)=BIG;
                                           
                                    Pnm(:,1)=xb.*temp5;
                                    temp5=temp5.*0+BIGI; 
                                                                      
                                    if m<nmax
                                       yb=xb;
                                       iyb=ixb;

                                       xb=(am(m+1).*tq).*yb;
                                       ixb=iyb;
                                       wb=abs(xb);

                                       wBIGSb=wb>=BIGS;
                                       wBIGSIb=wb<BIGSI;
                                       temp3(wBIGSb)=BIGI;
                                       temp3(wBIGSIb)=BIG;
                                       temp4(wBIGSb)=1;
                                       temp4(wBIGSIb)=-1;

                                       xb=xb.*temp3;
                                       ixb=ixb+temp4;
                                       temp3=temp2;
                                       temp4=temp4.*0;
                                       
                                       temp5(ixb==0)=1;
                                       temp5(ixb<-1)=0;
                                       %temp5(izb>=-1 & izb<0)=BIGI;
                                       %The condition "izb>=-1 & izb<0"
                                       %is useless, as "izb" is already
                                       %initialized as "izb=BIGI".
                                       temp5(ixb>0)=BIG;
                                           
                                       Pnm(:,2)=xb.*temp5;
                                       temp5=temp5.*0+BIGI; 

                                       for n=m+2:nmax
                                           idb=ixb-iyb;

                                           id0b=idb==0;
                                           id1b=idb==1;
                                           id_1b=idb==-1;
                                           idv1b=idb>1;
                                           
                                           temp1(id0b)=1;
                                           temp1(id1b)=1;
                                           temp2(id1b)=BIGI;
                                           temp1(id_1b)=BIGI;
                                           temp1(idv1b)=1;
                                           temp2(idv1b)=0;
                                           
                                           zzb=(am(n).*tq).*(xb.*temp1)-bm(n).*((yb.*q2).*temp2);
                                           izb=iyb;
                                           id0b_id1b_idv1b=id0b | id1b | idv1b;
                                           izb(id0b_id1b_idv1b)=ixb(id0b_id1b_idv1b);
                                           temp1=temp1.*0;
                                           temp2=temp1+1;

                                           wb=abs(zzb);

                                           wBIGSb=wb>=BIGS;
                                           wBIGSIb=wb<BIGSI;
                                           temp3(wBIGSb)=BIGI;
                                           temp3(wBIGSIb)=BIG;
                                           temp4(wBIGSb)=1;
                                           temp4(wBIGSIb)=-1;

                                           zzb=zzb.*temp3;
                                           izb=izb+temp4;
                                           temp3=temp2;
                                           temp4=temp4.*0;

                                           temp5(izb==0)=1;
                                           temp5(izb<-1)=0;
                                           %temp5(izb>=-1 & izb<0)=BIGI;
                                           %The condition "izb>=-1 & izb<0"
                                           %is useless, as "izb" is already
                                           %initialized as "izb=BIGI".
                                           temp5(izb>0)=BIG;
                                           
                                           Pnm(:,n-m+1)=zzb.*temp5;
                                           temp5=temp1+BIGI;   
                              
                                           yb=xb;
                                           iyb=ixb;
                                           xb=zzb;
                                           ixb=izb;
                                       end   
                                    end
                                end
                            end
                        end

                        %If nmin~=0
                        %======================================================
                        if nmin>nminGGM
                            end_idx=min([nmin-m nmin-nminGGM]);
                            if LNOFnmin==0
                                if por==1
                                    deltaCm(1:end_idx)=0;
                                end

                                if grav==1
                                    Cm(1:end_idx)=0;
                                end

                                if geoid==1
                                    HCm(1:end_idx)=0;
                                    HSm(1:end_idx)=0;
                                end

                                Sm(1:end_idx)=0;
                            elseif LNOFnmin==1
                                Pnm(:,1:end_idx)=0;
                            end
                        end
                        %======================================================

                        %% Computation of the first-order derivatives of modified fnALFs
                        if dALFs==1  
                            if volbaALFs==1 || volbaALFs==3
                                enm=sqrt((2*(m:nmax)+1).*((m:nmax)-m).*((m:nmax)+m)./(2*(m:nmax)-1));

                                if m==0 %Zonal modified dALFs
                                    dPnm(:,1)=0.*u;
                                    dPnm(:,2)=sqrt(3)*u.*q;
                                    dPnm(:,3:end)=-bsxfun(@times,2:nmax,tu).*Pnm(:,3:end)+(bsxfun(@times,enm(3:end),qu)).*Pnm(:,2:(end-1));
                                elseif m==nmax %Sectorial modified dALFs
                                    dPnm(:,1)=-m*(tu).*Pnm(:,1);
                                else
                                    dPnm(:,1)=-m*(tu).*Pnm(:,1); %Tesseral modified dALFs
                                    dPnm(:,2:(nmax-m+1))=-bsxfun(@times,(m+1):nmax,tu).*Pnm(:,2:(nmax-m+1))+bsxfun(@times,enm(2:end),qu).*Pnm(:,1:(nmax-m)); %Teseralne ALFs
                                end

                            elseif volbaALFs==2
                                enm=sqrt((2*(m:nmax)+1).*((m:nmax)-m).*((m:nmax)+m)./(2*(m:nmax)-1));

                                if m==0 %Zonal modified dALFs
                                    dPnm(:,1)=0.*u*1e-280;
                                    dPnm(:,2)=sqrt(3)*u.*q*1e-280;
                                    dPnm(:,3:end)=-bsxfun(@times,2:nmax,tu).*Pnm(:,3:end)+(bsxfun(@times,enm(3:end),qu)).*Pnm(:,2:(end-1));
                                elseif m==nmax
                                    dPnm(:,1)=-m*(tu).*Pnm(:,1); %Sectorial modified dALFs
                                else
                                    dPnm(:,1)=-m*(tu).*Pnm(:,1); %Tesseral modified dALFs
                                    dPnm(:,2:(nmax-m+1))=-bsxfun(@times,(m+1):nmax,tu).*Pnm(:,2:(nmax-m+1))+bsxfun(@times,enm(2:end),qu).*Pnm(:,1:(nmax-m)); %Teseralne dALFs
                                end
                            end                   

                            %Treatment of the dALFs singularity
                            dPnm(singdPnm,:)=0;
                            
                            if ddALFs==1 %If the second-order derivatives of the modified fnALFs are to be computed
                                if m==0 %Zonal modified ddALFs
                                    ddPnm=bsxfun(@times,tu,dPnm)-bsxfun(@times,(0:nmax).*((0:nmax)+1),Pnm);
                                else
                                    
                                    ddPnm(:,1:end-m)=bsxfun(@times,tu,dPnm(:,1:end-m))+bsxfun(@times,m^2./u.^2,Pnm(:,1:end-m))-bsxfun(@times,(m:nmax).*((m:nmax)+1),Pnm(:,1:end-m));
                                end
                                                                
                                %Treatment of the ddALFs singularity
                                ddPnm(singddPnm,:)=0;
                            end
                            
                        end                                     

                        cosla=cos(m*lambda);
                        sinla=sin(m*lambda);

                        if any(volbapar==8) || any(volbapar==9) || any(volbapar==14) || any(volbapar==15)
                            if m<2
                                coslaplus2=cos((m+2)*lambda);
                                sinlaplus2=sin((m+2)*lambda);
                                
                                coslaminus2=cos(0*lambda);
                                sinlaminus2=sin(0*lambda);
                            elseif m>nmax-2
                                coslaplus2=cos(0*lambda);
                                sinlaplus2=sin(0*lambda);
                                
                                coslaminus2=cos((m-2)*lambda);
                                sinlaminus2=sin((m-2)*lambda);
                            else
                                coslaplus2=cos((m+2)*lambda);
                                sinlaplus2=sin((m+2)*lambda);
                                
                                coslaminus2=cos((m-2)*lambda);
                                sinlaminus2=sin((m-2)*lambda);
                            end
                        end
                        
                        if any(volbapar==9) || any(volbapar==15)
                            if m<1
                                coslaplus1=cos((m+1)*lambda);
                                sinlaplus1=sin((m+1)*lambda);
                                
                                coslaminus1=cos(0*lambda);
                                sinlaminus1=sin(0*lambda);
                            elseif m>nmax-1
                                coslaplus1=cos(0*lambda);
                                sinlaplus1=sin(0*lambda);
                                
                                coslaminus1=cos((m-1)*lambda);
                                sinlaminus1=sin((m-1)*lambda);
                            else
                                coslaplus1=cos((m+1)*lambda);
                                sinlaplus1=sin((m+1)*lambda);
                                
                                coslaminus1=cos((m-1)*lambda);
                                sinlaminus1=sin((m-1)*lambda);
                            end
                        end
                        
                        m_max=max([m nminGGM]);
                        
                        %% Loop for 1:NF (number of computing functionals)                        
                        for i=1:pocetpar 

                            %Summation over n
                            if volbapar(i)==1         
                            elseif volbapar(i)==2 %Deflection of the vertical eta                                
                                 
                                for n=m_max:nmax
                                    eta=eta+(-deltaCm(n-m_max+1).*sinla+Sm(n-m_max+1).*cosla).*(m*Pnm(:,n-m+1));
                                end

                                if volbaALFs==2
                                    etaH=etaH.*u+eta;
                                    eta=0;
                                end

                            elseif volbapar(i)==3 %Deflection of the vertical xi

                                for n=m_max:nmax
                                    ksi=ksi+(deltaCm(n-m_max+1).*cosla+Sm(n-m_max+1).*sinla).*dPnm(:,n-m+1);
                                end
                                
                                if volbaALFs==2
                                    ksiH=ksiH.*u+ksi;
                                    ksi=0;
                                end

                            elseif volbapar(i)==4 %Deflection of the vertical Theta
                                
                                for n=m_max:nmax
                                    Teta=Teta+(-deltaCm(n-m_max+1).*sinla+Sm(n-m_max+1).*cosla).*(m*Pnm(:,n-m+1));
                                    Tksi=Tksi+(deltaCm(n-m_max+1).*cosla+Sm(n-m_max+1).*sinla).*dPnm(:,n-m+1);
                                end

                                if volbaALFs==2
                                    TetaH=TetaH.*u+Teta;
                                    Teta=0;
                                    TksiH=TksiH.*u+Tksi;
                                    Tksi=0;
                                end

                            elseif volbapar(i)==5 %Disturbing potential
                                
                                for n=m_max:nmax
                                    T=T+(deltaCm(n-m_max+1).*cosla+Sm(n-m_max+1).*sinla).*Pnm(:,n-m+1);
                                end

                                if volbaALFs==2
                                    TH=TH.*u+T;
                                    T=0;
                                end

                            elseif volbapar(i)==6 %Disturbing tensor Trr_Tpp_Tll
                             
                                for n=m_max:nmax
                                    deltaCmcosla=deltaCm(n-m_max+1).*cosla;
                                    Smsinla=Sm(n-m_max+1).*sinla;

                                    Trr=Trr+(n+1).*(n+2).*(deltaCmcosla+Smsinla).*Pnm(:,n-m+1);
                                    Tfifi=Tfifi+(deltaCmcosla+Smsinla).*ddPnm(:,n-m+1);
                                    Tll=Tll+(deltaCmcosla+Smsinla).*m.^2.*Pnm(:,n-m+1);
                                end

                                if volbaALFs==2
                                    TrrH=TrrH.*u+Trr;
                                    Trr=0;
                                    TfifiH=TfifiH.*u+Tfifi;
                                    Tfifi=0;
                                    TllH=TllH.*u+Tll;
                                    Tll=0;
                                end
                                
                            elseif volbapar(i)==7 %Disturbing tensor Trp_Trl_Tpl
                                                              
                                for n=m_max:nmax
                                    Smcosla=Sm(n-m_max+1).*cosla;
                                    deltaCmsinla=deltaCm(n-m_max+1).*sinla;

                                    Trfi=Trfi+(n+1).*(deltaCm(n-m_max+1).*cosla+Sm(n-m_max+1).*sinla).*dPnm(:,n-m+1);
                                    Trl=Trl+(n+1).*(Smcosla-deltaCmsinla).*(m*Pnm(:,n-m+1));
                                    Tfil=Tfil+(Smcosla-deltaCmsinla).*(m*dPnm(:,n-m+1));
                                end

                                if volbaALFs==2
                                    TrfiH=TrfiH.*u+Trfi;
                                    Trfi=0;
                                    TrlH=TrlH.*u+Trl;
                                    Trl=0;
                                    TfilH=TfilH.*u+Tfil;
                                    Tfil=0;
                                end
                                
                            elseif volbapar(i)==8 %Disturbing tensor Txx_Tyy_Tzz
                                
                                %Coefficients of the ALFs
                                if m==0
                                    anm=sqrt(2)/4.*sqrt((m:nmax).^2-(m+2-1).^2).*sqrt((m:nmax)+m+2).*sqrt((m:nmax)-(m+2)+2);
                                    bnm=((m:nmax)+m+1).*((m:nmax)+m+2)./2./(m+1);
                                    cnm=zeros(nmax+1,1);
                                elseif m==1
                                    anm=1/4.*sqrt((m:nmax).^2-(m+2-1).^2).*sqrt((m:nmax)+m+2).*sqrt((m:nmax)-(m+2)+2);
                                    bnm=((m:nmax)+m+1).*((m:nmax)+m+2)./2./(m+1);
                                    cnm=zeros(nmax,1);
                                elseif m==2
                                    anm=1/4.*sqrt((m:nmax).^2-(m+2-1).^2).*sqrt((m:nmax)+m+2).*sqrt((m:nmax)-(m+2)+2);
                                    bnm=((m:nmax).^2+m^2+3.*(m:nmax)+2)./2;
                                    cnm=sqrt(2)/4.*sqrt((m:nmax).^2-(m-2+1).^2).*sqrt((m:nmax)-(m-2)).*sqrt((m:nmax)+m-2+2);
                                else
                                    anm=1/4.*sqrt((m:nmax).^2-(m+2-1).^2).*sqrt((m:nmax)+m+2).*sqrt((m:nmax)-(m+2)+2);
                                    bnm=((m:nmax).^2+m^2+3.*(m:nmax)+2)./2;
                                    cnm=1/4.*sqrt((m:nmax).^2-(m-2+1).^2).*sqrt((m:nmax)-(m-2)).*sqrt((m:nmax)+m-2+2);
                                end
                                  
                                %Spherical harmonic coefficients of orders m+2 and m-2
                                if m==0 && nmax<=2
                                    deltaCplus2=deltaC(index((m+m_min-nminGGM):end)+m+2);
                                    Splus2=S(index((m+m_min-nminGGM):end)+m+2); 

                                    deltaCminus2=zeros(nmax-m+1,1);
                                    Sminus2=zeros(nmax-m+1,1);
                                elseif m<2 && nmax<=2
                                    deltaCplus2=zeros(nmax-m+1,1);
                                    Splus2=zeros(nmax-m+1,1);
                                    
                                    deltaCminus2=zeros(nmax-m+1,1);
                                    Sminus2=zeros(nmax-m+1,1);
                                elseif m<2
                                    deltaCplus2=deltaC(index((m+m_min-nminGGM):end)+m+2);
                                    Splus2=S(index((m+m_min-nminGGM):end)+m+2); 
                                        
                                    deltaCminus2=zeros(nmax-m+1,1);
                                    Sminus2=zeros(nmax-m+1,1); 
                                elseif m>nmax-2
                                    deltaCplus2=zeros(nmax-m+1,1);
                                    Splus2=zeros(nmax-m+1,1);
                                        
                                    deltaCminus2=deltaC(index((m+m_min-nminGGM):end)+m-2);
                                    Sminus2=S(index((m+m_min-nminGGM):end)+m-2); 
                                else
                                    deltaCplus2=deltaC(index((m+m_min-nminGGM):end)+m+2);
                                    Splus2=S(index((m+m_min-nminGGM):end)+m+2);
                                        
                                    deltaCminus2=deltaC(index((m+m_min-nminGGM):end)+m-2);
                                    Sminus2=S(index((m+m_min-nminGGM):end)+m-2); 
                                end
                                
                                %Summation over m                    
                                for n=m_max:nmax
                                    %Tzz
                                    Tzz=Tzz+(n+1).*(n+2).*(deltaCm(n-m_max+1).*cosla+Sm(n-m_max+1).*sinla).*Pnm(:,n-m+1);

                                    %Txx
                                    Txx1=(deltaCplus2(n-m_max+1).*coslaplus2+Splus2(n-m_max+1).*sinlaplus2).*anm(n-m+1);
                                    Txx2_vnutro=deltaCm(n-m_max+1).*cosla+Sm(n-m_max+1).*sinla;
                                    Txx2=Txx2_vnutro.*(bnm(n-m+1)-(n+1)*(n+2));
                                    Txx3=(deltaCminus2(n-m_max+1).*coslaminus2+Sminus2(n-m_max+1).*sinlaminus2).*cnm(n-m+1);
                                    Txx=Txx+(Txx1+Txx2+Txx3).*Pnm(:,n-m+1);

                                    %Tyy
                                    Tyy2=Txx2_vnutro.*bnm(n-m+1);
                                    Tyy=Tyy+(Txx1+Tyy2+Txx3).*Pnm(:,n-m+1); 
                                end

                                if volbaALFs==2
                                    TzzH=TzzH.*u+Tzz;
                                    Tzz=0;
                                    TxxH=TxxH.*u+Txx;
                                    Txx=0;
                                    TyyH=TyyH.*u+Tyy;
                                    Tyy=0;
                                end
                                
                            elseif volbapar(i)==9 %Disturbing tensor Txy_Txz_Tyz
                                                               
                                %Coefficients of the ALFs
                                if m==0
                                    dnm=1/4.*sqrt((2*(m:nmax)+1)./(2*(m:nmax)-1)).*sqrt(2).*sqrt((m:nmax).^2-(m+2-1).^2).*sqrt((m:nmax)+m+2).*sqrt((m:nmax)+m+2-2);
                                    gnm=zeros(nmax+1,1);
                                    hnm=gnm;
                                    
                                    betanm=((m:nmax)+2)./2.*sqrt(1+ones(1,nmax+1)).*sqrt((m:nmax)+m+1).*sqrt((m:nmax)-(m+1)+1);
                                    gamanm=gnm;
                                    
                                    minm=((m:nmax)+2)./2.*sqrt((2*(m:nmax)+1)./(2*(m:nmax)-1)).*sqrt(2).*sqrt((m:nmax)+m+1).*sqrt((m:nmax)+m+1-1);
                                    ninm=gnm;
                                elseif m==1
                                    dnm=1/4.*sqrt((2*(m:nmax)+1)./(2*(m:nmax)-1)).*sqrt((m:nmax).^2-(m+2-1).^2).*sqrt((m:nmax)+m+2).*sqrt((m:nmax)+m+2-2);
                                    gnm=-1/4*sqrt((2*(m:nmax)+1)./(2*(m:nmax)-1)).*sqrt((m:nmax)+1).*sqrt((m:nmax)-1).*((m:nmax)+2);
                                    hnm=zeros(nmax,1);
                                    
                                    betanm=((m:nmax)+2)./2.*sqrt((m:nmax)+m+1).*sqrt((m:nmax)-(m+1)+1);
                                    gamanm=-((m:nmax)+2).*sqrt((m:nmax).*((m:nmax)+1)./2);
                                    
                                    minm=((m:nmax)+2)./2.*sqrt((2*(m:nmax)+1)./(2*(m:nmax)-1)).*sqrt((m:nmax)+m+1).*sqrt((m:nmax)+m+1-1);
                                    ninm=gnm;
                                elseif m==2
                                    dnm=1/4.*sqrt((2*(m:nmax)+1)./(2*(m:nmax)-1)).*sqrt((m:nmax).^2-(m+2-1).^2).*sqrt((m:nmax)+m+2).*sqrt((m:nmax)+m+2-2);
                                    gnm=-m/2*sqrt((2*(m:nmax)+1)./(2*(m:nmax)-1)).*sqrt((m:nmax)+m).*sqrt((m:nmax)-m);
                                    hnm=zeros(nmax-1,1);
                                    
                                    betanm=((m:nmax)+2)./2.*sqrt((m:nmax)+m+1).*sqrt((m:nmax)-(m+1)+1);
                                    gamanm=-((m:nmax)+2)./2.*sqrt((m:nmax)-(m-1)).*sqrt((m:nmax)+m-1+1);
                                    
                                    minm=((m:nmax)+2)./2.*sqrt((2*(m:nmax)+1)./(2*(m:nmax)-1)).*sqrt((m:nmax)+m+1).*sqrt((m:nmax)+m+1-1);
                                    ninm=((m:nmax)+2)./2.*sqrt((2*(m:nmax)+1)./(2*(m:nmax)-1)).*sqrt((m:nmax)-(m-1)).*sqrt((m:nmax)-(m-1)-1);
                                elseif m==3
                                    dnm=1/4.*sqrt((2*(m:nmax)+1)./(2*(m:nmax)-1)).*sqrt((m:nmax).^2-(m+2-1).^2).*sqrt((m:nmax)+m+2).*sqrt((m:nmax)+m+2-2);
                                    gnm=-m/2*sqrt((2*(m:nmax)+1)./(2*(m:nmax)-1)).*sqrt((m:nmax)+m).*sqrt((m:nmax)-m);
                                    hnm=-1/4*sqrt((2*(m:nmax)+1)./(2*(m:nmax)-1)).*sqrt((m:nmax)-3).*sqrt((m:nmax)-2).*sqrt((m:nmax)-1).*sqrt((m:nmax)+2);
                                    
                                    betanm=((m:nmax)+2)./2.*sqrt((m:nmax)+m+1).*sqrt((m:nmax)-(m+1)+1);
                                    gamanm=-((m:nmax)+2)./2.*sqrt((m:nmax)-(m-1)).*sqrt((m:nmax)+m-1+1);
                                    
                                    minm=((m:nmax)+2)./2.*sqrt((2*(m:nmax)+1)./(2*(m:nmax)-1)).*sqrt((m:nmax)+m+1).*sqrt((m:nmax)+m+1-1);
                                    ninm=((m:nmax)+2)./2.*sqrt((2*(m:nmax)+1)./(2*(m:nmax)-1)).*sqrt((m:nmax)-(m-1)).*sqrt((m:nmax)-(m-1)-1);
                                else
                                    dnm=1/4.*sqrt((2*(m:nmax)+1)./(2*(m:nmax)-1)).*sqrt((m:nmax).^2-(m+2-1).^2).*sqrt((m:nmax)+m+2).*sqrt((m:nmax)+m+2-2);
                                    gnm=-m/2*sqrt((2*(m:nmax)+1)./(2*(m:nmax)-1)).*sqrt((m:nmax)+m).*sqrt((m:nmax)-m);
                                    hnm=-1/4*sqrt((2*(m:nmax)+1)./(2*(m:nmax)-1)).*sqrt((m:nmax).^2-(m-2+1).^2).*sqrt((m:nmax)-(m-2)).*sqrt((m:nmax)-(m-2)-2);
                                    
                                    betanm=((m:nmax)+2)./2.*sqrt((m:nmax)+m+1).*sqrt((m:nmax)-(m+1)+1);
                                    gamanm=-((m:nmax)+2)./2.*sqrt((m:nmax)-(m-1)).*sqrt((m:nmax)+m-1+1);
                                    
                                    minm=((m:nmax)+2)./2.*sqrt((2*(m:nmax)+1)./(2*(m:nmax)-1)).*sqrt((m:nmax)+m+1).*sqrt((m:nmax)+m+1-1);
                                    ninm=((m:nmax)+2)./2.*sqrt((2*(m:nmax)+1)./(2*(m:nmax)-1)).*sqrt((m:nmax)-(m-1)).*sqrt((m:nmax)-(m-1)-1);
                                end
                                
                                %Spherical harmonic coefficients of orders m+2, m-2, m+1 and m-1
                                if m==0 && nmax<=2
                                    deltaCplus2=deltaC(index((m+m_min-nminGGM):end)+m+2);
                                    Splus2=S(index((m+m_min-nminGGM):end)+m+2);
                                    
                                    deltaCminus2=zeros(nmax-m+1,1);
                                    Sminus2=zeros(nmax-m+1,1);
                                    
                                    deltaCplus1=deltaC(index((m+m_min-nminGGM):end)+m+1);
                                    Splus1=S(index((m+m_min-nminGGM):end)+m+1);
                                    
                                    deltaCminus1=zeros(nmax-m+1,1);
                                    Sminus1=zeros(nmax-m+1,1);
                                elseif m<2 && nmax<=2
                                    deltaCplus2=zeros(nmax-m+1,1);
                                    Splus2=zeros(nmax-m+1,1);
                                    
                                    deltaCminus2=zeros(nmax-m+1,1);
                                    Sminus2=zeros(nmax-m+1,1);
                                    
                                    deltaCplus1=deltaC(index((m+m_min-nminGGM):end)+m+1);
                                    Splus1=S(index((m+m_min-nminGGM):end)+m+1);
                                    
                                    deltaCminus1=deltaC(index((m+m_min-nminGGM):end)+m-1);
                                    Sminus1=S(index((m+m_min-nminGGM):end)+m-1);
                                elseif m<2                                   
                                    deltaCplus2=deltaC(index((m+m_min-nminGGM):end)+m+2);
                                    Splus2=S(index((m+m_min-nminGGM):end)+m+2); 
                                        
                                    deltaCminus2=zeros(nmax-m+1,1);
                                    Sminus2=deltaCminus2; 
                                        
                                    deltaCplus1=deltaC(index((m+m_min-nminGGM):end)+m+1);
                                    Splus1=S(index((m+m_min-nminGGM):end)+m+1);
                                        
                                    if m==0
                                        deltaCminus1=zeros(nmax+1,1);
                                        Sminus1=deltaCminus1;
                                    else
                                        deltaCminus1=deltaC(index((m+m_min-nminGGM):end)+m-1);
                                        Sminus1=S(index((m+m_min-nminGGM):end)+m-1);
                                    end
                                elseif m>nmax-2
                                    deltaCplus2=zeros(nmax-m+1,1);
                                    Splus2=deltaCplus2;
                                        
                                    deltaCminus2=deltaC(index((m+m_min-nminGGM):end)+m-2);
                                    Sminus2=S(index((m+m_min-nminGGM):end)+m-2); 
                                        
                                    if m==nmax
                                        deltaCplus1=zeros(nmax-m+1,1);
                                        Splus1=deltaCplus1;
                                    else
                                        deltaCplus1=deltaC(index((m+m_min-nminGGM):end)+m+1);
                                        Splus1=S(index((m+m_min-nminGGM):end)+m+1);
                                    end
                                        
                                    deltaCminus1=deltaC(index((m+m_min-nminGGM):end)+m-1);
                                    Sminus1=S(index((m+m_min-nminGGM):end)+m-1); 
                                else
                                    deltaCplus2=deltaC(index((m+m_min-nminGGM):end)+m+2);
                                    Splus2=S(index((m+m_min-nminGGM):end)+m+2);
                                        
                                    deltaCminus2=deltaC(index((m+m_min-nminGGM):end)+m-2);
                                    Sminus2=S(index((m+m_min-nminGGM):end)+m-2); 
                                        
                                    deltaCplus1=deltaC(index((m+m_min-nminGGM):end)+m+1);
                                    Splus1=S(index((m+m_min-nminGGM):end)+m+1);
                                        
                                    deltaCminus1=deltaC(index((m+m_min-nminGGM):end)+m-1);
                                    Sminus1=S(index((m+m_min-nminGGM):end)+m-1); 
                                end
                                
                                PTnm=[zeros(length_fi,1) Pnm];
                                
                                %Summation over n                          
                                for n=m_max:nmax
                                    %Txy
                                    Txy1=(-deltaCplus2(n-m_max+1).*sinlaplus2+Splus2(n-m_max+1).*coslaplus2).*dnm(n-m+1);
                                    Txy2=(-deltaCm(n-m_max+1).*sinla+Sm(n-m_max+1).*cosla).*gnm(n-m+1);
                                    Txy3=(-deltaCminus2(n-m_max+1).*sinlaminus2+Sminus2(n-m_max+1).*coslaminus2).*hnm(n-m+1);
                                    Txy=Txy+(Txy1+Txy2+Txy3).*PTnm(:,n-m+1).*q;

                                    %Txz
                                    Txz1=(deltaCplus1(n-m_max+1).*coslaplus1+Splus1(n-m_max+1).*sinlaplus1).*betanm(n-m+1);
                                    Txz2=(deltaCminus1(n-m_max+1).*coslaminus1+Sminus1(n-m_max+1).*sinlaminus1).*gamanm(n-m+1);
                                    Txz=Txz+(Txz1+Txz2).*Pnm(:,n-m+1);

                                    %Tyz
                                    Tyz1=(-deltaCplus1(n-m_max+1).*sinlaplus1+Splus1(n-m_max+1).*coslaplus1).*minm(n-m+1);
                                    Tyz2=(-deltaCminus1(n-m_max+1).*sinlaminus1+Sminus1(n-m_max+1).*coslaminus1).*ninm(n-m+1);
                                    Tyz=Tyz+(Tyz1+Tyz2).*PTnm(:,n-m+1).*q;  
                                end

                                if volbaALFs==2
                                    TxyH=TxyH.*u+Txy;
                                    Txy=0;
                                    TxzH=TxzH.*u+Txz;
                                    Txz=0;
                                    TyzH=TyzH.*u+Tyz;
                                    Tyz=0;
                                end
                                
                            elseif volbapar(i)==10 %Geoid undulation
                                   
                                for n=m_max:nmax
                                    N1c=N1c+(deltaCm(n-m_max+1).*cosla+Sm(n-m_max+1).*sinla).*Pnm(:,n-m+1);
                                        
                                    % When computing H, there is no
                                    % dumping factor (R./r).^n,
                                    % therefore the matrix Pnm has to be
                                    % devided by 1./((R./r).^n), since
                                    % Pnm is the matrix of the MODIFIED
                                    % fnALFS
                                    H=H+1./(R./r).^n.*(HCm(n-m_max+1).*cosla+HSm(n-m_max+1).*sinla).*Pnm(:,n-m+1); 
                                end
                                    
                                if volbaALFs==2
                                    N1cH=N1cH.*u+N1c;
                                    N1c=0;
                                    HH=HH.*u+H;
                                    H=0;
                                end                            
                                                        
                            elseif volbapar(i)==11 %Gravitational potential
                               
                                for n=m_max:nmax 
                                    V=V+(Cm(n-m_max+1).*cosla+Sm(n-m_max+1).*sinla).*Pnm(:,n-m+1);
                                end
                                
                                if volbaALFs==2    
                                    VH=VH.*u+V;
                                    V=0;
                                end

                            elseif volbapar(i)==12 %Gravitational tensor Vrr_Vpp_Vll
                                
                                for n=m_max:nmax
                                    Cmcosla=Cm(n-m_max+1).*cosla;
                                    Smsinla=Sm(n-m_max+1).*sinla;

                                    Vrr=Vrr+(n+1).*(n+2).*(Cmcosla+Smsinla).*Pnm(:,n-m+1);
                                    Vfifi=Vfifi+(Cmcosla+Smsinla).*ddPnm(:,n-m+1);
                                    Vll=Vll+(Cmcosla+Smsinla).*m.^2.*Pnm(:,n-m+1);
                                end

                                if volbaALFs==2
                                    VrrH=VrrH.*u+Vrr;
                                    Vrr=0;
                                    VfifiH=VfifiH.*u+Vfifi;
                                    Vfifi=0;
                                    VllH=VllH.*u+Vll;
                                    Vll=0;
                                end
                                
                            elseif volbapar(i)==13 %Gravitational tensor Vrp_Vrl_Vpl
                                
                                for n=m_max:nmax
                                    Smcosla=Sm(n-m_max+1).*cosla;
                                    Cmsinla=Cm(n-m_max+1).*sinla;

                                    Vrfi=Vrfi+(n+1).*(Cm(n-m_max+1).*cosla+Sm(n-m_max+1).*sinla).*dPnm(:,n-m+1);
                                    Vrl=Vrl+(n+1).*(Smcosla-Cmsinla).*(m*Pnm(:,n-m+1));
                                    Vfil=Vfil+(Smcosla-Cmsinla).*(m*dPnm(:,n-m+1));
                                end

                                if volbaALFs==2
                                    VrfiH=VrfiH.*u+Vrfi;
                                    Vrfi=0;
                                    VrlH=VrlH.*u+Vrl;
                                    Vrl=0;
                                    VfilH=VfilH.*u+Vfil;
                                    Vfil=0;
                                end
                                
                            elseif volbapar(i)==14 %Gravitational tensor Vxx_Vyy_Vzz
                                
                                %Coefficients of the ALFs
                                if m==0
                                    anm=sqrt(2)/4.*sqrt((m:nmax).^2-(m+2-1).^2).*sqrt((m:nmax)+m+2).*sqrt((m:nmax)-(m+2)+2);
                                    bnm=((m:nmax)+m+1).*((m:nmax)+m+2)./2./(m+1);
                                    cnm=zeros(nmax+1,1);
                                elseif m==1
                                    anm=1/4.*sqrt((m:nmax).^2-(m+2-1).^2).*sqrt((m:nmax)+m+2).*sqrt((m:nmax)-(m+2)+2);
                                    bnm=((m:nmax)+m+1).*((m:nmax)+m+2)./2./(m+1);
                                    cnm=zeros(nmax,1);
                                elseif m==2
                                    anm=1/4.*sqrt((m:nmax).^2-(m+2-1).^2).*sqrt((m:nmax)+m+2).*sqrt((m:nmax)-(m+2)+2);
                                    bnm=((m:nmax).^2+m^2+3.*(m:nmax)+2)./2;
                                    cnm=sqrt(2)/4.*sqrt((m:nmax).^2-(m-2+1).^2).*sqrt((m:nmax)-(m-2)).*sqrt((m:nmax)+m-2+2);
                                else
                                    anm=1/4.*sqrt((m:nmax).^2-(m+2-1).^2).*sqrt((m:nmax)+m+2).*sqrt((m:nmax)-(m+2)+2);
                                    bnm=((m:nmax).^2+m^2+3.*(m:nmax)+2)./2;
                                    cnm=1/4.*sqrt((m:nmax).^2-(m-2+1).^2).*sqrt((m:nmax)-(m-2)).*sqrt((m:nmax)+m-2+2);
                                end
                                  
                                %Spherical harmonic coefficients of orders m+2 a m-2
                                if m==0 && nmax<=2
                                    Cplus2=C(index((m+m_min-nminGGM):end)+m+2);
                                    Splus2=S(index((m+m_min-nminGGM):end)+m+2); 

                                    Cminus2=zeros(nmax-m+1,1);
                                    Sminus2=zeros(nmax-m+1,1);
                                elseif m<2 && nmax<=2
                                    Cplus2=zeros(nmax-m+1,1);
                                    Splus2=zeros(nmax-m+1,1);
                                    
                                    Cminus2=zeros(nmax-m+1,1);
                                    Sminus2=zeros(nmax-m+1,1);
                                elseif m<2
                                    Cplus2=C(index((m+m_min-nminGGM):end)+m+2);
                                    Splus2=S(index((m+m_min-nminGGM):end)+m+2); 
                                        
                                    Cminus2=zeros(nmax-m+1,1);
                                    Sminus2=zeros(nmax-m+1,1); 
                                elseif m>nmax-2
                                    Cplus2=zeros(nmax-m+1,1);
                                    Splus2=zeros(nmax-m+1,1);
                                        
                                    Cminus2=C(index((m+m_min-nminGGM):end)+m-2);
                                    Sminus2=S(index((m+m_min-nminGGM):end)+m-2); 
                                else
                                    Cplus2=C(index((m+m_min-nminGGM):end)+m+2);
                                    Splus2=S(index((m+m_min-nminGGM):end)+m+2);
                                        
                                    Cminus2=C(index((m+m_min-nminGGM):end)+m-2);
                                    Sminus2=S(index((m+m_min-nminGGM):end)+m-2); 
                                end
                                
                                %Summation over n                            
                                for n=m_max:nmax
                                    %Vzz
                                    Vzz=Vzz+(n+1).*(n+2).*(Cm(n-m_max+1).*cosla+Sm(n-m_max+1).*sinla).*Pnm(:,n-m+1);

                                    %Vxx
                                    Vxx1=(Cplus2(n-m_max+1).*coslaplus2+Splus2(n-m_max+1).*sinlaplus2).*anm(n-m+1);
                                    Vxx2_vnutro=Cm(n-m_max+1).*cosla+Sm(n-m_max+1).*sinla;
                                    Vxx2=Vxx2_vnutro.*(bnm(n-m+1)-(n+1)*(n+2));
                                    Vxx3=(Cminus2(n-m_max+1).*coslaminus2+Sminus2(n-m_max+1).*sinlaminus2).*cnm(n-m+1);
                                    Vxx=Vxx+(Vxx1+Vxx2+Vxx3).*Pnm(:,n-m+1);

                                    %Vyy
                                    Vyy2=Vxx2_vnutro.*bnm(n-m+1);
                                    Vyy=Vyy+(Vxx1+Vyy2+Vxx3).*Pnm(:,n-m+1); 
                                end

                                if volbaALFs==2
                                    VzzH=VzzH.*u+Vzz;
                                    Vzz=0;
                                    VxxH=VxxH.*u+Vxx;
                                    Vxx=0;
                                    VyyH=VyyH.*u+Vyy;
                                    Vyy=0;
                                end
                                
                            elseif volbapar(i)==15 %Gravitational tensor Vxy_Vxz_Vyz
                                
                                %Coefficietns of the ALFs
                                if m==0
                                    dnm=1/4.*sqrt((2*(m:nmax)+1)./(2*(m:nmax)-1)).*sqrt(2).*sqrt((m:nmax).^2-(m+2-1).^2).*sqrt((m:nmax)+m+2).*sqrt((m:nmax)+m+2-2);
                                    gnm=zeros(nmax+1,1);
                                    hnm=gnm;
                                    
                                    betanm=((m:nmax)+2)./2.*sqrt(1+ones(1,nmax+1)).*sqrt((m:nmax)+m+1).*sqrt((m:nmax)-(m+1)+1);
                                    gamanm=gnm;
                                    
                                    minm=((m:nmax)+2)./2.*sqrt((2*(m:nmax)+1)./(2*(m:nmax)-1)).*sqrt(2).*sqrt((m:nmax)+m+1).*sqrt((m:nmax)+m+1-1);
                                    ninm=gnm;
                                elseif m==1
                                    dnm=1/4.*sqrt((2*(m:nmax)+1)./(2*(m:nmax)-1)).*sqrt((m:nmax).^2-(m+2-1).^2).*sqrt((m:nmax)+m+2).*sqrt((m:nmax)+m+2-2);
                                    gnm=-1/4*sqrt((2*(m:nmax)+1)./(2*(m:nmax)-1)).*sqrt((m:nmax)+1).*sqrt((m:nmax)-1).*((m:nmax)+2);
                                    hnm=zeros(nmax,1);
                                    
                                    betanm=((m:nmax)+2)./2.*sqrt((m:nmax)+m+1).*sqrt((m:nmax)-(m+1)+1);
                                    gamanm=-((m:nmax)+2).*sqrt((m:nmax).*((m:nmax)+1)./2);
                                    
                                    minm=((m:nmax)+2)./2.*sqrt((2*(m:nmax)+1)./(2*(m:nmax)-1)).*sqrt((m:nmax)+m+1).*sqrt((m:nmax)+m+1-1);
                                    ninm=gnm;
                                elseif m==2
                                    dnm=1/4.*sqrt((2*(m:nmax)+1)./(2*(m:nmax)-1)).*sqrt((m:nmax).^2-(m+2-1).^2).*sqrt((m:nmax)+m+2).*sqrt((m:nmax)+m+2-2);
                                    gnm=-m/2*sqrt((2*(m:nmax)+1)./(2*(m:nmax)-1)).*sqrt((m:nmax)+m).*sqrt((m:nmax)-m);
                                    hnm=zeros(nmax-1,1);
                                    
                                    betanm=((m:nmax)+2)./2.*sqrt((m:nmax)+m+1).*sqrt((m:nmax)-(m+1)+1);
                                    gamanm=-((m:nmax)+2)./2.*sqrt((m:nmax)-(m-1)).*sqrt((m:nmax)+m-1+1);
                                    
                                    minm=((m:nmax)+2)./2.*sqrt((2*(m:nmax)+1)./(2*(m:nmax)-1)).*sqrt((m:nmax)+m+1).*sqrt((m:nmax)+m+1-1);
                                    ninm=((m:nmax)+2)./2.*sqrt((2*(m:nmax)+1)./(2*(m:nmax)-1)).*sqrt((m:nmax)-(m-1)).*sqrt((m:nmax)-(m-1)-1);
                                elseif m==3
                                    dnm=1/4.*sqrt((2*(m:nmax)+1)./(2*(m:nmax)-1)).*sqrt((m:nmax).^2-(m+2-1).^2).*sqrt((m:nmax)+m+2).*sqrt((m:nmax)+m+2-2);
                                    gnm=-m/2*sqrt((2*(m:nmax)+1)./(2*(m:nmax)-1)).*sqrt((m:nmax)+m).*sqrt((m:nmax)-m);
                                    hnm=-1/4*sqrt((2*(m:nmax)+1)./(2*(m:nmax)-1)).*sqrt((m:nmax)-3).*sqrt((m:nmax)-2).*sqrt((m:nmax)-1).*sqrt((m:nmax)+2);
                                    
                                    betanm=((m:nmax)+2)./2.*sqrt((m:nmax)+m+1).*sqrt((m:nmax)-(m+1)+1);
                                    gamanm=-((m:nmax)+2)./2.*sqrt((m:nmax)-(m-1)).*sqrt((m:nmax)+m-1+1);
                                    
                                    minm=((m:nmax)+2)./2.*sqrt((2*(m:nmax)+1)./(2*(m:nmax)-1)).*sqrt((m:nmax)+m+1).*sqrt((m:nmax)+m+1-1);
                                    ninm=((m:nmax)+2)./2.*sqrt((2*(m:nmax)+1)./(2*(m:nmax)-1)).*sqrt((m:nmax)-(m-1)).*sqrt((m:nmax)-(m-1)-1);
                                else
                                    dnm=1/4.*sqrt((2*(m:nmax)+1)./(2*(m:nmax)-1)).*sqrt((m:nmax).^2-(m+2-1).^2).*sqrt((m:nmax)+m+2).*sqrt((m:nmax)+m+2-2);
                                    gnm=-m/2*sqrt((2*(m:nmax)+1)./(2*(m:nmax)-1)).*sqrt((m:nmax)+m).*sqrt((m:nmax)-m);
                                    hnm=-1/4*sqrt((2*(m:nmax)+1)./(2*(m:nmax)-1)).*sqrt((m:nmax).^2-(m-2+1).^2).*sqrt((m:nmax)-(m-2)).*sqrt((m:nmax)-(m-2)-2);
                                    
                                    betanm=((m:nmax)+2)./2.*sqrt((m:nmax)+m+1).*sqrt((m:nmax)-(m+1)+1);
                                    gamanm=-((m:nmax)+2)./2.*sqrt((m:nmax)-(m-1)).*sqrt((m:nmax)+m-1+1);
                                    
                                    minm=((m:nmax)+2)./2.*sqrt((2*(m:nmax)+1)./(2*(m:nmax)-1)).*sqrt((m:nmax)+m+1).*sqrt((m:nmax)+m+1-1);
                                    ninm=((m:nmax)+2)./2.*sqrt((2*(m:nmax)+1)./(2*(m:nmax)-1)).*sqrt((m:nmax)-(m-1)).*sqrt((m:nmax)-(m-1)-1);
                                end
                                
                                %Spherical harmonic coefficients of orders m+2, m-2, m+1 and m-1
                                if m==0 && nmax<=2
                                    Cplus2=C(index((m+m_min-nminGGM):end)+m+2);
                                    Splus2=S(index((m+m_min-nminGGM):end)+m+2);
                                    
                                    Cminus2=zeros(nmax-m+1,1);
                                    Sminus2=zeros(nmax-m+1,1);
                                    
                                    Cplus1=C(index((m+m_min-nminGGM):end)+m+1);
                                    Splus1=S(index((m+m_min-nminGGM):end)+m+1);
                                    
                                    Cminus1=zeros(nmax-m+1,1);
                                    Sminus1=zeros(nmax-m+1,1);
                                elseif m<2 && nmax<=2
                                    Cplus2=zeros(nmax-m+1,1);
                                    Splus2=zeros(nmax-m+1,1);
                                    
                                    Cminus2=zeros(nmax-m+1,1);
                                    Sminus2=zeros(nmax-m+1,1);
                                    
                                    Cplus1=C(index((m+m_min-nminGGM):end)+m+1);
                                    Splus1=S(index((m+m_min-nminGGM):end)+m+1);
                                    
                                    Cminus1=C(index((m+m_min-nminGGM):end)+m-1);
                                    Sminus1=S(index((m+m_min-nminGGM):end)+m-1);
                                elseif m<2
                                    Cplus2=C(index((m+m_min-nminGGM):end)+m+2);
                                    Splus2=S(index((m+m_min-nminGGM):end)+m+2); 
                                        
                                    Cminus2=zeros(nmax-m+1,1);
                                    Sminus2=Cminus2; 
                                        
                                    Cplus1=C(index((m+m_min-nminGGM):end)+m+1);
                                    Splus1=S(index((m+m_min-nminGGM):end)+m+1);
                                       
                                    if m==0
                                        Cminus1=zeros(nmax+1,1);
                                        Sminus1=Cminus1;
                                    else
                                        Cminus1=C(index((m+m_min-nminGGM):end)+m-1);
                                        Sminus1=S(index((m+m_min-nminGGM):end)+m-1);
                                    end
                                elseif m>nmax-2
                                    Cplus2=zeros(nmax-m+1,1);
                                    Splus2=Cplus2;
                                        
                                    Cminus2=C(index((m+m_min-nminGGM):end)+m-2);
                                    Sminus2=S(index((m+m_min-nminGGM):end)+m-2); 
                                        
                                    if m==nmax
                                        Cplus1=zeros(nmax-m+1,1);
                                        Splus1=Cplus1;
                                    else
                                        Cplus1=C(index((m+m_min-nminGGM):end)+m+1);
                                        Splus1=S(index((m+m_min-nminGGM):end)+m+1);
                                    end
                                        
                                    Cminus1=C(index((m+m_min-nminGGM):end)+m-1);
                                    Sminus1=S(index((m+m_min-nminGGM):end)+m-1); 
                                else
                                    Cplus2=C(index((m+m_min-nminGGM):end)+m+2);
                                    Splus2=S(index((m+m_min-nminGGM):end)+m+2);
                                        
                                    Cminus2=C(index((m+m_min-nminGGM):end)+m-2);
                                    Sminus2=S(index((m+m_min-nminGGM):end)+m-2); 
                                       
                                    Cplus1=C(index((m+m_min-nminGGM):end)+m+1);
                                    Splus1=S(index((m+m_min-nminGGM):end)+m+1);
                                        
                                    Cminus1=C(index((m+m_min-nminGGM):end)+m-1);
                                    Sminus1=S(index((m+m_min-nminGGM):end)+m-1); 
                                end
                                
                                PTnm=[zeros(length_fi,1) Pnm];
                                
                                %Summation over n                            
                                for n=m_max:nmax
                                    %Vxy
                                    Vxy1=(-Cplus2(n-m_max+1).*sinlaplus2+Splus2(n-m_max+1).*coslaplus2).*dnm(n-m+1);
                                    Vxy2=(-Cm(n-m_max+1).*sinla+Sm(n-m_max+1).*cosla).*gnm(n-m+1);
                                    Vxy3=(-Cminus2(n-m_max+1).*sinlaminus2+Sminus2(n-m_max+1).*coslaminus2).*hnm(n-m+1);
                                    Vxy=Vxy+(Vxy1+Vxy2+Vxy3).*PTnm(:,n-m+1).*q;

                                    %Vxz
                                    Vxz1=(Cplus1(n-m_max+1).*coslaplus1+Splus1(n-m_max+1).*sinlaplus1).*betanm(n-m+1);
                                    Vxz2=(Cminus1(n-m_max+1).*coslaminus1+Sminus1(n-m_max+1).*sinlaminus1).*gamanm(n-m+1);
                                    Vxz=Vxz+(Vxz1+Vxz2).*Pnm(:,n-m+1);

                                    %Vyz
                                    Vyz1=(-Cplus1(n-m_max+1).*sinlaplus1+Splus1(n-m_max+1).*coslaplus1).*minm(n-m+1);
                                    Vyz2=(-Cminus1(n-m_max+1).*sinlaminus1+Sminus1(n-m_max+1).*coslaminus1).*ninm(n-m+1);
                                    Vyz=Vyz+(Vyz1+Vyz2).*PTnm(:,n-m+1).*q;  
                                end

                                if volbaALFs==2
                                    VxyH=VxyH.*u+Vxy;
                                    Vxy=0;
                                    VxzH=VxzH.*u+Vxz;
                                    Vxz=0;
                                    VyzH=VyzH.*u+Vyz;
                                    Vyz=0;
                                end
                                
                            elseif volbapar(i)==16 %Gravity vector gX_gY_gZ
                               
                                for n=m_max:nmax 
                                    Cmcosla=Cm(n-m_max+1).*cosla;
                                    Smsinla=Sm(n-m_max+1).*sinla;
                                        
                                    Wr=Wr+(n+1).*(Cmcosla+Smsinla).*Pnm(:,n-m+1);
                                    Wlambda=Wlambda+(-Cm(n-m_max+1).*sinla+Sm(n-m_max+1).*cosla).*(m*Pnm(:,n-m+1)); 
                                    Wfi=Wfi+(Cmcosla+Smsinla).*dPnm(:,n-m+1); 
                                end

                                if volbaALFs==2
                                    WrH=WrH.*u+Wr;
                                    Wr=0;
                                    WlambdaH=WlambdaH.*u+Wlambda;
                                    Wlambda=0;
                                    WfiH=WfiH.*u+Wfi;
                                    Wfi=0;
                                end

                            elseif volbapar(i)==17 %Gravity sa
                            
                                for n=m_max:nmax  
                                    g_sa=g_sa+(n+1).*(Cm(n-m_max+1).*cosla+Sm(n-m_max+1).*sinla).*Pnm(:,n-m+1);
                                end

                                if volbaALFs==2
                                    g_saH=g_saH.*u+g_sa;
                                    g_sa=0;
                                end

                            elseif volbapar(i)==18 %Gravity potential
                                                               
                                for n=m_max:nmax 
                                    W=W+(Cm(n-m_max+1).*cosla+Sm(n-m_max+1).*sinla).*Pnm(:,n-m+1);
                                end

                                if volbaALFs==2
                                    WH=WH.*u+W;
                                    W=0;
                                end

                            elseif volbapar(i)==19 %Gravity anomaly sa  
                                                               
                                for n=m_max:nmax 
                                    anomalia_sa=anomalia_sa+(n-1).*(deltaCm(n-m_max+1).*cosla+Sm(n-m_max+1).*sinla).*Pnm(:,n-m+1);
                                end

                                if volbaALFs==2
                                    anomalia_saH=anomalia_saH.*u+anomalia_sa;
                                    anomalia_sa=0;
                                end

                            elseif volbapar(i)==20 %Gravity disturbance
                              
                                for n=m_max:nmax 
                                    Cmcosla=Cm(n-m_max+1).*cosla;
                                    Smsinla=Sm(n-m_max+1).*sinla;
                                        
                                    Wr=Wr+(n+1).*(Cmcosla+Smsinla).*Pnm(:,n-m+1);
                                    Wlambda=Wlambda+(-Cm(n-m_max+1).*sinla+Sm(n-m_max+1).*cosla).*(m*Pnm(:,n-m+1)); 
                                    Wfi=Wfi+(Cmcosla+Smsinla).*dPnm(:,n-m+1);  
                                end

                                if m==0
                                    for n=m_max:nmax
                                        Ur=Ur+(n+1).*CElm(n-m_max+1).*cosla.*Pnm(:,n-m+1); 
                                        Ufi=Ufi+CElm(n-m_max+1).*cosla.*dPnm(:,n-m+1); 
                                    end
                                end

                                if volbaALFs==2
                                    WrH=WrH.*u+Wr;
                                    Wr=0;
                                    WlambdaH=WlambdaH.*u+Wlambda;
                                    Wlambda=0;
                                    WfiH=WfiH.*u+Wfi;
                                    Wfi=0;
                                    UrH=UrH.*u+Ur;
                                    Ur=0;
                                    UfiH=UfiH.*u+Ufi;
                                    Ufi=0;
                                end

                            elseif volbapar(i)==21 %Gravity disturbance sa

                                for n=m_max:nmax 
                                    porucha_sa=porucha_sa+(n+1).*(deltaCm(n-m_max+1).*cosla+Sm(n-m_max+1).*sinla).*Pnm(:,n-m+1);
                                end
                                
                                if volbaALFs==2
                                    porucha_saH=porucha_saH.*u+porucha_sa;
                                    porucha_sa=0;
                                end

                            elseif volbapar(i)==22 %Height anomaly Ell
                                                             
                                for n=m_max:nmax 
                                    zetaEl=zetaEl+(deltaCm(n-m_max+1).*cosla+Sm(n-m_max+1).*sinla).*Pnm(:,n-m+1);
                                end

                                if volbaALFs==2
                                    zetaElH=zetaElH.*u+zetaEl;
                                    zetaEl=0;
                                end
                                
                            elseif volbapar(i)==23 %Height anomaly
                                                               
                                for n=m_max:nmax 
                                    deltaCmcosla_Smsinla_Pnm=(deltaCm(n-m_max+1).*cosla+Sm(n-m_max+1).*sinla).*Pnm(:,n-m+1);
                                    zeta_N1c=zeta_N1c+deltaCmcosla_Smsinla_Pnm;
                                        
                                    % When computing H, there is no
                                    % dumping factor (R./r).^n,
                                    % therefore the matrix Pnm has to be
                                    % devided by 1./((R./r).^n), since
                                    % Pnm is the matrix of the MODIFIED
                                    % fnALFS
                                    zeta_H=zeta_H+1./(R./r).^n.*(HCm(n-m_max+1).*cosla+HSm(n-m_max+1).*sinla).*Pnm(:,n-m+1); 
                                    zeta_dg=zeta_dg+(n+1).*deltaCmcosla_Smsinla_Pnm;
                                end
                                    
                                if volbaALFs==2
                                    zeta_N1cH=zeta_N1cH.*u+zeta_N1c;
                                    zeta_N1c=0;
                                    zeta_HH=zeta_HH.*u+zeta_H;
                                    zeta_H=0;
                                    zeta_dgH=zeta_dgH.*u+zeta_dg;
                                    zeta_dg=0;
                                end
                                
                            elseif volbapar(i)==24 %Second radial derivative of disturbing potential

                                for n=m_max:nmax    
                                    T_rr=T_rr+(n+1).*(n+2).*(deltaCm(n-m_max+1).*cosla+Sm(n-m_max+1).*sinla).*Pnm(:,n-m+1);
                                end

                                if volbaALFs==2
                                    T_rrH=T_rrH.*u+T_rr;
                                    T_rr=0;
                                end

                            elseif volbapar(i)==25 %Second radial derivative of gravity potential

                               for n=m_max:nmax 
                                    Wrr=Wrr+(n+1).*(n+2).*(Cm(n-m_max+1).*cosla+Sm(n-m_max+1).*sinla).*Pnm(:,n-m+1);
                               end

                               if volbaALFs==2
                                   WrrH=WrrH.*u+Wrr;
                                   Wrr=0;
                               end
                                    
                            end                                              
                        end


                    end
                    
                    %Update of the progress bar
                    if GUI==1 %If working with the GUI
                        set(progressbar,'string',...
                            'Progress: Matrix multiplications...',...
                            'fontsize',8); drawnow;
                    end
                    
                    if GUI==0 && Status_bar==1 %If working without the GUI
                        fprintf('\n')
                        fprintf('Progress: Matrix multiplications...\n')
                    end

                    clear Pnm dPnm ddPnm Cm Sm C CElm deltaC deltaCm ...
                        S u t q q2 index sinla cosla ...
                        deltaCmcosla_Smsinla_Pnm Cmcosla deltaCmcosla ...
                        Smcosla deltaCmsinla Smsinla enm singdPnm ...
                        qu tu singddPnm 

                    %Computation of the normal gravity for eta, xi, Theta,
                    %N, zeta el, zeta
                    if any(volbapar==2) || any(volbapar==3) || any(volbapar==4) || any(volbapar==10) || any(volbapar==22) || any(volbapar==23)                        
                        bEl=aEl*sqrt(1-eEl^2);
                        EEl=sqrt(aEl^2-bEl^2);
                        
                        %Computation of ellipsoidal harmonic coordinates
                        [X,Y,Z]=ell2cart(fi,lambda,h,[aEl eEl]);
                        ugama=sqrt((X.*X+Y.*Y+Z.*Z-EEl.^2).*(1./2+1./2.*sqrt(1+((4.*EEl.^2.*Z.*Z)./(X.*X+Y.*Y+Z.*Z-EEl.^2).^2))));                        
                        betagama=atan((Z.*sqrt(ugama.^2+EEl.^2))./(ugama.*sqrt(X.*X+Y.*Y)));
                        clear X Y Z
                        
                        wgama=sqrt((ugama.^2+EEl^2*sin(betagama).^2)./(ugama.^2+EEl^2));
                        qgama=1/2*((1+(3*ugama.^2)./EEl^2).*atan(EEl./ugama)-3*ugama./EEl);
                        qgama_=3*(1+(ugama.^2)./EEl^2).*(1-ugama./EEl.*atan(EEl./ugama))-1;
                        qgama0=1/2*((1+(3*bEl^2)/EEl^2)*atan(EEl/bEl)-3*bEl/EEl);
                        
                        gamau=-1./wgama.*(GMEl./(ugama.^2+EEl.^2)+(omegaEl.^2.*aEl.^2.*EEl)./(ugama.^2+EEl.^2).*(qgama_./qgama0).*(1./2.*sin(betagama).^2-1./6)-omegaEl.^2.*ugama.*cos(betagama).^2);
                        gamabeta=-1./wgama.*(-(omegaEl.^2.*aEl.^2.*qgama)./(sqrt(ugama.^2+EEl.^2).*qgama0)+omegaEl.^2.*sqrt(ugama.^2+EEl.^2)).*sin(betagama).*cos(betagama);
                        
                        clear ugama betagama wgama qgama qgama_ qgama0
                        
                        gamaP=sqrt(gamau.^2+gamabeta.^2);
                        
                        clear gamau gamabeta
                    end
                    
                    %% Final computation of functionals of the geopotential 
                    for i=1:pocetpar 
                        if volbapar(i)==1                
                        elseif volbapar(i)==2 %Deflection of the vertical eta

                            if volbaALFs==1 || volbaALFs==3
                                Pg=-(GM./(r.^2.*gamaP.*cos(fiG)).*(sum(eta,2)))*(180/pi)*3600;
                            elseif volbaALFs==2
                                Pg=-(GM./(r.^2.*gamaP.*cos(fiG)).*(etaH*1e280))*(180/pi)*3600;
                            end

                            Pg(fi==pi/2 | fi==-pi/2)=0;
                            
                            clear eta etaH

                        elseif volbapar(i)==3 %Deflection of the vertical xi

                            if volbaALFs==1 || volbaALFs==3
                                Pg=-(GM./(r.^2.*gamaP).*(sum(ksi,2)))*(180/pi)*3600;
                            elseif volbaALFs==2
                                Pg=-(GM./(r.^2.*gamaP).*(ksiH*1e280))*(180/pi)*3600;
                            end

                            clear ksi ksiH

                        elseif volbapar(i)==4 %Deflection of the vertical Theta

                            if volbaALFs==1 || volbaALFs==3
                                Teta=-(GM./(r.^2.*gamaP.*cos(fiG)).*(sum(Teta,2)))*(180/pi)*3600;
                                Tksi=-(GM./(r.^2.*gamaP).*(sum(Tksi,2)))*(180/pi)*3600;
                                
                                Talfa=atan2(Teta,Tksi);
                                Talfa(Talfa<0)=Talfa(Talfa<0)+2*pi;
                            
                                Pg=[sqrt(Teta.^2+Tksi.^2) 180/pi*(Talfa)];
                            
                                clear Teta Tksi Talfa
                                
                                Pg(fi==pi/2 | fi==-pi/2,:)=0;
                                
                            elseif volbaALFs==2
                                TetaH=-(GM./(r.^2.*gamaP.*cos(fiG)).*(TetaH*1e280))*(180/pi)*3600;
                                TksiH=-(GM./(r.^2.*gamaP).*(TksiH*1e280))*(180/pi)*3600;
                                
                                Talfa=atan2(TetaH,TksiH);
                                Talfa(Talfa<0)=Talfa(Talfa<0)+2*pi;

                                Pg=[sqrt(TetaH.^2+TksiH.^2) 180/pi*(Talfa)];

                                clear TetaH TksiH Talfa
                                
                                Pg(fi==pi/2 | fi==-pi/2,:)=0;
                            end

                        elseif volbapar(i)==5 %Disturbing potential

                            if volbaALFs==1 || volbaALFs==3
                                Pg=(GM./r.*(sum(T,2)));
                            elseif volbaALFs==2
                                Pg=(GM./r.*(TH*1e280));
                            end

                            clear T TH

                        elseif volbapar(i)==6 %Disturbing tensor Trr_Tpp_Tll
                            
                            if volbaALFs==1 || volbaALFs==3
                                Trr=(GM./r.^3.*(sum(Trr,2)))*10^9;
                                Tfifi=(GM./r.^3.*(sum(Tfifi,2)))*10^9;
                                Tll=(GM./(r.^3.*cos(fiG).^2).*(sum(Tll,2)))*10^9;
                                
                                Pg=Trr;
                                clear Trr
                                Pg=[Pg Tfifi];
                                clear Tfifi
                                Pg=[Pg -Tll];
                                clear Tll
                            elseif volbaALFs==2
                                TrrH=(GM./r.^3.*(TrrH*1e280))*10^9;
                                TfifiH=(GM./r.^3.*(TfifiH*1e280))*10^9;
                                TllH=(GM./(r.^3.*cos(fiG).^2).*(TllH*1e280))*10^9;
                                
                                Pg=TrrH;
                                clear TrrH
                                Pg=[Pg TfifiH];
                                clear TfifiH
                                Pg=[Pg -TllH];
                                clear TllH
                            end
                            
                        elseif volbapar(i)==7 %Disturbing tensor Trp_Trl_Tpl
                            
                            if volbaALFs==1 || volbaALFs==3
                                Trfi=(GM./r.^3.*(sum(Trfi,2)))*10^9;
                                Trl=(GM./(r.^3.*cos(fiG)).*(sum(Trl,2)))*10^9;
                                Tfil=(GM./(r.^3.*cos(fiG)).*(sum(Tfil,2)))*10^9;
                                
                                Pg=-Trfi;
                                clear Trfi
                                Pg=[Pg -Trl];
                                clear Trl
                                Pg=[Pg Tfil];
                                clear Tfil
                            elseif volbaALFs==2
                                TrfiH=(GM./r.^3.*(TrfiH*1e280))*10^9;
                                TrlH=(GM./(r.^3.*cos(fiG)).*(TrlH*1e280))*10^9;
                                TfilH=(GM./(r.^3.*cos(fiG)).*(TfilH*1e280))*10^9;
                                
                                Pg=-TrfiH;
                                clear TrfiH
                                Pg=[Pg -TrlH];
                                clear TrlH
                                Pg=[Pg TfilH];
                                clear TfilH
                            end
                            
                        elseif volbapar(i)==8 %Disturbing tensor Txx_Tyy_Tzz
                            
                            clear deltaCplus2 Splus2 coslaplus2 sinlaplus2 ...
                                deltaCminus2 Sminus2 coslaminus2 sinlaminus2 ...
                                Txx1 Txx2 Txx3 Tyy2 Txx2_vnutro anm bnm cnm
                            
                            if volbaALFs==1 || volbaALFs==3
                                Tzz=(GM./r.^3.*(sum(Tzz,2)))*10^9;
                                Txx=(GM./r.^3.*(sum(Txx,2)))*10^9;
                                Tyy=(GM./r.^3.*(sum(Tyy,2)))*10^9;
                                
                                Pg=[Txx]; %#ok<*NBRAK>
                                clear Txx
                                Pg=[Pg -Tyy];
                                clear Tyy
                                Pg=[Pg Tzz];
                                clear Tzz
                            elseif volbaALFs==2
                                TzzH=(GM./r.^3.*(TzzH*1e280))*10^9;
                                TxxH=(GM./r.^3.*(TxxH*1e280))*10^9;
                                TyyH=(GM./r.^3.*(TyyH*1e280))*10^9;
                                
                                Pg=[TxxH];
                                clear TxxH
                                Pg=[Pg -TyyH];
                                clear TyyH
                                Pg=[Pg TzzH];
                                clear TzzH
                            end
                            
                        elseif volbapar(i)==9 %Disturbing tensor Txy_Txz_Tyz
                            
                            clear deltaCplus2 Splus2 coslaplus2 sinlaplus2 ...
                                deltaCminus2 Sminus2 coslaminus2 sinlaminus2 ...
                                deltaCplus1 Splus1 coslaplus1 sinlaplus1 ...
                                deltaCminus1 Sminus1 coslaminus1 sinlaminus1 ...
                                PTnm Txy1 Txy2 Txy3 Txz1 Txz2 Tyz1 Tyz2 ...
                                dnm gnm hnm betanm gamanm minm ninm
                            
                            if volbaALFs==1 || volbaALFs==3
                                Txy=(GM./r.^3.*(sum(Txy,2)))*10^9;
                                Txz=(GM./r.^3.*(sum(Txz,2)))*10^9;
                                Tyz=(GM./r.^3.*(sum(Tyz,2)))*10^9;
                                
                                Pg=[Txy];
                                clear Txy
                                Pg=[Pg Txz];
                                clear Txz
                                Pg=[Pg Tyz];
                                clear Tyz
                            elseif volbaALFs==2
                                TxyH=(GM./r.^3.*(TxyH*1e280))*10^9;
                                TxzH=(GM./r.^3.*(TxzH*1e280))*10^9;
                                TyzH=(GM./r.^3.*(TyzH*1e280))*10^9;
                                
                                Pg=[TxyH];
                                clear TxyH
                                Pg=[Pg TxzH];
                                clear TxzH
                                Pg=[Pg TyzH];
                                clear TyzH
                            end

                        elseif volbapar(i)==10 %Geoid undulation
                            
                            clear HC HCm HS HSm
                            
                            G=6.67259*10^-11; %Newtonian gravitational constant (Moritz, 2000, Geodetic reference system 1980)
                            ro=2670; %Density of the crust

                            if volbaALFs==1 || volbaALFs==3
                                N1c=(GM./r.*(sum(N1c,2)))./gamaP;
                                H(H<0)=H(H<0)*0; %H is set to zero in the areas of oceans and seas
                                Pg=N1c-(2*pi*G*ro*H.^2)./gamaP;
                                
                                clear N1c H
                                
                            elseif volbaALFs==2
                                N1cH=(GM./r.*(N1cH*1e280))./gamaP;
                                HH=HH*1e280;
                                HH(HH<0)=HH(HH<0)*0;
                                Pg=N1cH-(2*pi*G*ro*HH.^2)./gamaP;
                                
                                clear N1cH HH
                            end
                             
                        elseif volbapar(i)==11 %Gravitational potential

                            if volbaALFs==1 || volbaALFs==3
                                Pg=(GM./r.*(sum(V,2)));
                            elseif volbaALFs==2
                                Pg=(GM./r.*(VH*1e280));
                            end

                            clear V VH

                        elseif volbapar(i)==12 %Gravitational tensor Vrr_Vpp_Vll
                            
                            if volbaALFs==1 || volbaALFs==3
                                Vrr=(GM./r.^3.*(sum(Vrr,2)))*10^9;
                                Vfifi=(GM./r.^3.*(sum(Vfifi,2)))*10^9;
                                Vll=(GM./(r.^3.*cos(fiG).^2).*(sum(Vll,2)))*10^9;
                                
                                Pg=Vrr;
                                clear Vrr
                                Pg=[Pg Vfifi];
                                clear Vfifi
                                Pg=[Pg -Vll];
                                clear Vll
                            elseif volbaALFs==2
                                VrrH=(GM./r.^3.*(VrrH*1e280))*10^9;
                                VfifiH=(GM./r.^3.*(VfifiH*1e280))*10^9;
                                VllH=(GM./(r.^3.*cos(fiG).^2).*(VllH*1e280))*10^9;
                                
                                Pg=VrrH;
                                clear VrrH
                                Pg=[Pg VfifiH];
                                clear VfifiH
                                Pg=[Pg -VllH];
                                clear VllH
                            end
                            
                        elseif volbapar(i)==13 %Gravitational tensor Vrp_Vrl_Vpl
                            
                            if volbaALFs==1 || volbaALFs==3
                                Vrfi=(GM./r.^3.*(sum(Vrfi,2)))*10^9;
                                Vrl=(GM./(r.^3.*cos(fiG)).*(sum(Vrl,2)))*10^9;
                                Vfil=(GM./(r.^3.*cos(fiG)).*(sum(Vfil,2)))*10^9;
                                
                                Pg=-Vrfi;
                                clear Vrfi
                                Pg=[Pg -Vrl];
                                clear Vrl
                                Pg=[Pg Vfil];
                                clear Vfil
                            elseif volbaALFs==2
                                VrfiH=(GM./r.^3.*(VrfiH*1e280))*10^9;
                                VrlH=(GM./(r.^3.*cos(fiG)).*(VrlH*1e280))*10^9;
                                VfilH=(GM./(r.^3.*cos(fiG)).*(VfilH*1e280))*10^9;
                                
                                Pg=-VrfiH;
                                clear VrfiH
                                Pg=[Pg -VrlH];
                                clear VrlH
                                Pg=[Pg VfilH];
                                clear VfilH
                            end
                            
                        elseif volbapar(i)==14 %Gravitational tensor Vxx_Vyy_Vzz
                            
                            clear Cplus2 Splus2 coslaplus2 sinlaplus2 ...
                                Cminus2 Sminus2 coslaminus2 sinlaminus2 ...
                                Vxx1 Vxx2 Vxx3 Vyy2 Vxx2_vnutro anm bnm cnm
                            
                            if volbaALFs==1 || volbaALFs==3
                                Vzz=(GM./r.^3.*(sum(Vzz,2)))*10^9;
                                Vxx=(GM./r.^3.*(sum(Vxx,2)))*10^9;
                                Vyy=(GM./r.^3.*(sum(Vyy,2)))*10^9;
                                
                                Pg=[Vxx];
                                clear Vxx
                                Pg=[Pg -Vyy];
                                clear Vyy
                                Pg=[Pg Vzz];
                                clear Vzz
                            elseif volbaALFs==2
                                VzzH=(GM./r.^3.*(VzzH*1e280))*10^9;
                                VxxH=(GM./r.^3.*(VxxH*1e280))*10^9;
                                VyyH=(GM./r.^3.*(VyyH*1e280))*10^9;
                                
                                Pg=[VxxH];
                                clear VxxH
                                Pg=[Pg -VyyH];
                                clear VyyH
                                Pg=[Pg VzzH];
                                clear VzzH
                            end
                            
                        elseif volbapar(i)==15 %%Gravitational tensor Vxy_Vxz_Vyz
                            
                            clear Cplus2 Splus2 coslaplus2 sinlaplus2 ...
                                Cminus2 Sminus2 coslaminus2 sinlaminus2 ...
                                Cplus1 Splus1 coslaplus1 sinlaplus1 ...
                                Cminus1 Sminus1 coslaminus1 sinlaminus1 ...
                                PTnm Vxy1 Vxy2 Vxy3 Vxz1 Vxz2 Vyz1 Vyz2 ...
                                dnm gnm hnm betanm gamanm minm ninm
                            
                            if volbaALFs==1 || volbaALFs==3
                                Vxy=(GM./r.^3.*(sum(Vxy,2)))*10^9;
                                Vxz=(GM./r.^3.*(sum(Vxz,2)))*10^9;
                                Vyz=(GM./r.^3.*(sum(Vyz,2)))*10^9;
                                
                                Pg=[Vxy];
                                clear Vxy
                                Pg=[Pg Vxz];
                                clear Vxz
                                Pg=[Pg Vyz];
                                clear Vyz
                            elseif volbaALFs==2
                                VxyH=(GM./r.^3.*(VxyH*1e280))*10^9;
                                VxzH=(GM./r.^3.*(VxzH*1e280))*10^9;
                                VyzH=(GM./r.^3.*(VyzH*1e280))*10^9;
                                
                                Pg=[VxyH];
                                clear VxyH
                                Pg=[Pg VxzH];
                                clear VxzH
                                Pg=[Pg VyzH];
                                clear VyzH
                            end
                            
                        elseif volbapar(i)==16 %Gravity vector gX_gY_gZ

                            if volbaALFs==1 || volbaALFs==3
                                Wr=(-GM./r.^2.*(sum(Wr,2))+omegaEl^2.*r.*(cos(fiG).^2));
                                Wlambda=-GM./r.*Wlambda;
                                Wfi=(GM./r.*(sum(Wfi,2))-omegaEl^2.*(r.^2).*cos(fiG).*sin(fiG));
                                                                
                                Pg=(1./r.*Wfi)*10^5;
                                clear Wfi
                                Pg=[Pg (Wlambda./(r.*cos(fiG)))*10^5];                                
                                clear Wlambda
                                Pg=[Pg Wr*10^5];
                                clear Wr
                            elseif volbaALFs==2
                                WrH=(-GM./r.^2.*(WrH*1e280)+omegaEl^2.*r.*(cos(fiG).^2));
                                WlambdaH=-GM./r.*WlambdaH*1e280;
                                WfiH=(GM./r.*WfiH*1e280-omegaEl^2.*(r.^2).*cos(fiG).*sin(fiG));
                                
                                Pg=(1./r.*WfiH)*10^5;
                                clear WfiH
                                Pg=[Pg (WlambdaH./(r.*cos(fiG)))*10^5];
                                clear WlambdaH
                                Pg=[Pg WrH*10^5];
                                clear WrH
                            end

                        elseif volbapar(i)==17 %Gravity sa

                            if volbaALFs==1 || volbaALFs==3
                                g_sa=(-GM./r.^2.*(sum(g_sa,2))+omegaEl^2.*r.*(cos(fiG).^2));
                                Pg=sqrt(g_sa.^2)*10^5;
                                
                                clear g_sa
                            elseif volbaALFs==2
                                g_saH=(-GM./r.^2.*(g_saH*1e280)+omegaEl^2.*r.*(cos(fiG).^2));
                                Pg=sqrt(g_saH.^2)*10^5;
                                
                                clear g_saH
                            end

                        elseif volbapar(i)==18 %Gravity potential

                            if volbaALFs==1 || volbaALFs==3
                                Pg=(GM./r.*(sum(W,2)))+1/2*omegaEl.^2.*r.^2.*cos(fiG).^2;
                            elseif volbaALFs==2
                                Pg=(GM./r.*(WH*1e280))+1/2*omegaEl.^2.*r.^2.*cos(fiG).^2;
                            end

                            clear W WH

                        elseif volbapar(i)==19 %Gravity anomaly sa

                            if volbaALFs==1 || volbaALFs==3
                                Pg=(GM./r.^2.*(sum(anomalia_sa,2)))*10^5;
                            elseif volbaALFs==2
                                Pg=(GM./r.^2.*(anomalia_saH*1e280))*10^5;
                            end
                            
                            clear anomalia_sa anomalia_saH
                            
                        elseif volbapar(i)==20 %Gravity disturbance

                            if volbaALFs==1 || volbaALFs==3
                                Wr=-GM./r.^2.*sum(Wr,2)+omegaEl^2.*r.*(cos(fiG).^2);
                                Wlambda=GM./r.*Wlambda;
                                Wfi=(GM./r.*sum(Wfi,2)-omegaEl^2.*(r.^2).*cos(fiG).*sin(fiG));
                                Ur=-GM./r.^2.*sum(Ur,2)+omegaEl^2.*r.*(cos(fiG).^2);
                                Ufi=(GM./r.*(sum(Ufi,2))-omegaEl^2.*(r.^2).*cos(fiG).*sin(fiG));
                                
                                g=sqrt(Wr.^2+(Wlambda./(r.*cos(fiG))).^2+(1./r.*Wfi).^2);
                                clear Wr Wlambda Wfi
                                
                                gama_GGM=sqrt(Ur.^2+(1./r.*Ufi).^2);
                                clear Ur Ufi

                                Pg=(g-gama_GGM)*10^5;
                                clear g gama_GGM                                
                            elseif volbaALFs==2
                                WrH=-GM./r.^2.*WrH*1e280+omegaEl^2.*r.*(cos(fiG).^2);
                                WlambdaH=GM./r.*WlambdaH*1e280;
                                WfiH=(GM./r.*WfiH*1e280-omegaEl^2.*(r.^2).*cos(fiG).*sin(fiG));
                                UrH=-GM./r.^2.*UrH*1e280+omegaEl^2.*r.*(cos(fiG).^2);
                                UfiH=(GM./r.*UfiH*1e280-omegaEl^2.*(r.^2).*cos(fiG).*sin(fiG));
                                
                                g=sqrt(WrH.^2+(WlambdaH./(r.*cos(fiG))).^2+(1./r.*WfiH).^2);
                                clear WrH WlambdaH WfiH
                                
                                gama_GGM=sqrt(UrH.^2+(1./r.*UfiH).^2);
                                clear UrH UfiH

                                Pg=(g-gama_GGM)*10^5;                                
                                clear g gama_GGM
                            end
                            
                        elseif volbapar(i)==21 %Gravity disturbance sa

                            if volbaALFs==1 || volbaALFs==3
                                Pg=(GM./r.^2.*(sum(porucha_sa,2)))*10^5;
                            elseif volbaALFs==2
                                Pg=(GM./r.^2.*(porucha_saH*1e280))*10^5;
                            end

                            clear porucha_sa porucha_saH
                        elseif volbapar(i)==22 %Height anomaly Ell  

                            if volbaALFs==1 || volbaALFs==3
                                Pg=(GM./(r.*gamaP).*(sum(zetaEl,2)));
                            elseif volbaALFs==2
                                Pg=(GM./(r.*gamaP).*(zetaElH*1e280));
                            end
                            
                            clear zetaEl zetaElH

                        elseif volbapar(i)==23 %Height anomaly
                            
                            clear HC HCm HS HSm
                            
                            G=6.67259*10^-11; %Newtonian gravitational constant (Moritz, 2000, Geodetic reference system 1980)
                            ro=2670; %Density of the crust

                            if volbaALFs==1 || volbaALFs==3
                                zeta_zetaEl=(GM./r.*(sum(zeta_N1c,2)))./gamaP;
                                zeta_N1c=(GM./r.*(sum(zeta_N1c,2)))./gamaP;
                                zeta_dg=(GM./r.^2.*(sum(zeta_dg,2)));
                                
                                zeta_H(zeta_H<0)=zeta_H(zeta_H<0)*0; %H is set to zero in the areas of oceans and seas                           
                                zeta_N=zeta_N1c-(2*pi*G*ro*zeta_H.^2)./gamaP;

                                Pg=zeta_zetaEl-zeta_dg.*((zeta_H+zeta_N)./gamaP);
                            
                                clear zeta_zetaEl zeta_dg zeta_H zeta_N zeta_N1c
                            elseif volbaALFs==2
                                zeta_zetaElH=(GM./r.*(zeta_N1cH*1e280))./gamaP;
                                zeta_N1cH=(GM./r.*(zeta_N1cH*1e280))./gamaP;
                                zeta_dgH=(GM./r.^2.*(zeta_dgH*1e280));
                                zeta_HH=zeta_HH*1e280;
                                
                                zeta_HH(zeta_HH<0)=zeta_HH(zeta_HH<0)*0;                          
                                zeta_NH=zeta_N1cH-(2*pi*G*ro*zeta_HH.^2)./gamaP;

                                Pg=zeta_zetaElH-zeta_dgH.*((zeta_HH+zeta_NH)./gamaP);
                            
                                clear zeta_zetaElH zeta_dgH zeta_HH zeta_NH zeta_N1cH
                            end
                            
                        elseif volbapar(i)==24 %Second radial derivative of disturbing potential

                            if volbaALFs==1 || volbaALFs==3
                                Pg=(GM./r.^3.*(sum(T_rr,2)))*10^9;
                            elseif volbaALFs==2
                                Pg=(GM./r.^3.*(T_rrH*1e280))*10^9;
                            end
                            
                            clear T_rr T_rrH

                        elseif volbapar(i)==25 %Second radial derivative of gravity potential

                            if volbaALFs==1 || volbaALFs==3
                                Pg=(GM./r.^3.*(sum(Wrr,2))+omegaEl^2.*cos(fiG).^2)*10^9;
                            elseif volbaALFs==2
                                Pg=(GM./r.^3.*(WrrH*1e280)+omegaEl^2.*cos(fiG).^2)*10^9;
                            end

                            clear Wrr WrrH
                        end

                        if i==1
                            P=Pg;
                            clear Pg
                        else
                            P=[P Pg];
                            if i==pocetpar
                                clear Pg
                            end
                        end
                    end 
                    
                    clear q q2 gamaP
                    
                    %Update of the progress bar
                    set(progressbar,'string','','fontsize',8); drawnow;
                    
                end
            end

            %% Computation of commission error            
            if STD==1
                
                tic %Start clock to measure computation time
                
                if GUI==0
                    if strcmp(nmax,'nmaxGGM')
                        use_nmax_GGM=1;
                    else
                        use_nmax_GGM=0;
                    end
                end
                
                %Warning for the computation of modified fnALFs
                if volbaALFs==2 || volbaALFs==3
                    volbaALFs=1;
                    if GUI==1 %If working with the GUI
                        warn1=warndlg('To compute commission error of functionals of the geopotential, only the standard forward column method can be used for computing fnALFs. After clicking OK, this method will be applied in the following computations.');
                        waitfor(warn1);
                    else
                        warning('To compute commission error of functionals of the geopotential, only the standard forward column method can be used for computing fnALFs. This method will be applied in the following computations.'); %#ok<WNTAG>
                    end
                end
                
                %Check if commission error of tensors is to be computed in
                %the LNOF (not allowed)
                if any(volbapar==8) || any(volbapar==9) || any(volbapar==14) || any(volbapar==15)
                    if GUI==1 %If working with the GUI
                        errordlg('Commission error of gravitational or disturbing tensor in the LNOF cannot be computed.',...
                            'Error in calculated parameters and output selection');
                    end
                    error('Commission error of gravitational or disturbing tensor in the LNOF cannot be computed.')
                end               

                %Error, if the input file has not been input
                if isempty(GGMcovname)
                    if GUI==1 %If working with the GUI
                        errordlg('Please input a file with the variance-covariance matrix of a geopotential model.',...
                            'Error in point type selection');
                    end
                    error('Please input a file with the variance-covariance matrix of a geopotential model.')
                end

                %Update of the progress bar
                if GUI==1 %If working with the GUI
                    set(findobj('tag','hlasky'),'string',...
                            'Loading covariance matrix...','fontsize',8,...
                            'foregroundcolor','k'); drawnow;
                end
                
                if GUI==0 && Status_bar==1 %If working without the GUI
                    fprintf('Loading covariance matrix...\n')
                end

                %Loading variance-covariance matrix
                if ~exist([GGMcovadresar,GGMcovname],'file') %Check whether the input file with the covariance matrix exists
                    if GUI==1 %Working with the GUI
                        errordlg('The entered file with the variance-covariance matrix does not exist.',...
                            'Error in point type selection');
                    end
                    error('The entered file with the variance-covariance matrix does not exist.')
                end
                if strcmp(GGMcovname(find(GGMcovname=='.'):end),'.mat') %If binary MAT file has been imported
                    covmat=importdata([GGMcovadresar,GGMcovname]);
                    nmaxGGMcov=max(covmat(:,2)); %Maximum degree of variance-covariance matrix
          
                    CSnm=covmat(:,1:3);
                    covmat=covmat(:,4:end);
                else %If ASCII file has been imported
                    fid=fopen([GGMcovadresar,GGMcovname]); 
                    covmat=fscanf(fid,'%f',[1 inf]);
                    fclose(fid);         

                    %Transformation of variance-covariance matrix from
                    %raw vector into the matrix
                    covmat=[0 0 0 0 0 0 covmat];

                    rows_GGM=length(covmat);
                    for i=0:2147483640
                        x=i^2-sum(1:(i-1));
                        if x==rows_GGM
                            rozmer=i;
                            break    	
                        end
                    end

                    korene=roots([1 2 -rozmer]);
                    nmaxGGMcov=korene(korene>0);

                    TM=triu(ones(nmaxGGMcov^2+2*nmaxGGMcov),0); %Lower triangular matrix with ones and zeros
                    TM(TM==1)=covmat;
                    covmat=TM'; %Lower triangular matrix with variances and covariances (raws 1:3 and columns 1:3 contain zeros)
                    
                    clear TM
                    
                    covmat=covmat(4:end,:);
  
                    if covmat(1,3)==0 && covmat(2,3)==1 %Sorted primarily according to degrees
                        if GUI==1 %If working with the GUI
                            errordlg('Wrong format of the input variance-covariance matrix.',...
                                'Geopotential model and reference system selection');
                        end
                        error('Wrong format of the input variance-covariance matrix.')
                    end
                                        
                    CSnm=covmat(:,1:3);
                    covmat=covmat(:,4:end);
                end                
                   
                [rows_covmat_or,cols_covmat_or]=size(covmat);
                if rows_covmat_or~=cols_covmat_or
                    if GUI==1 %If working with the GUI
                        errordlg('Wrong format of the input variance-covariance matrix.',...
                            'Geopotential model and reference system selection');
                    end
                    error('Wrong format of the input variance-covariance matrix.')
                end
                
                %Value of nmin and error messages
                if nmin<2
                    if GUI==1 %If working with the GUI
                        errordlg('To compute commission error of a functional, the value of "nmin" must be at least 2.',...
                            'Error in geopotential model and reference system selection')
                    end
                    error('To compute commission error of a functional, the value of "nmin" must be at least 2.')
                elseif nmin>nmaxGGMcov
                    if GUI==1 %If working with the GUI
                        errordlg('The value of "nmin" exceedes nmax value of the variance-covariance matrix of the GGM.',....
                            'Error in geopotential model and reference system selection')
                    end
                    error('The value of "nmin" exceedes nmax value of the variance-covariance matrix of the GGM.')
                end
                if isnan(nmin)==1
                    if GUI==1 %If working with the GUI
                        errordlg('Please input the "nmin" value.',...
                            'Error in geopotential model and reference system selection')
                    end
                    error('Please input the "nmin" value.')
                end

                nmin0=0; %Logical 0
                nmin1=0; %Logical 0
                
                %Value of nmax and error messages 
                if use_nmax_GGM==1 %If working with the GUI and maximum value of the GGM is used automatically
                    nmax=nmaxGGMcov;
                end
                if nmax>nmaxGGMcov
                    if GUI==1 %If working with the GUI
                        errordlg('The entered value of "nmax" exceedes "nmax" of the variance-covariance matrix.',...
                            'Error in geopotential model and reference system selection')
                    end
                    error('The entered value of "nmax" exceedes "nmax" of the variance-covariance matrix.')
                elseif nmax<2
                    if GUI==1 %If working with the GUI
                        errordlg('The value of "nmax" must be at least 2.',...
                            'Error in geopotential model and reference system selection')
                    end
                    error('The value of "nmax" must be at least 2.')
                elseif nmin>nmax
                    if GUI==1 %If working with the GUI
                        errordlg('The value of "nmin" cannot be larger than nmax value.',...
                            'Error in geopotential model and reference system selection')
                    end
                    error('The value of "nmin" cannot be larger than nmax value.')
                end
                if isnan(nmax)==1
                    if GUI==1 %If working with the GUI
                        errordlg('Please input the nmax value.',...
                            'Error in geopotential model and reference system selection')
                    end
                    error('Please input the nmax value.')
                end
               
                %Deleting of variances and covariaces if nmax < nmaxGGM
                if nmax<nmaxGGMcov
                    zmazat=CSnm(:,2)>nmax;
                    covmat(zmazat,:)=[];
                    covmat(:,zmazat)=[];
                end

                clear CSnm zmazat

                [rows_cov,cols_cov]=size(covmat);
                covmat=covmat'.*triu(ones(rows_cov),1)+covmat; %Creating full variances-covariance matrix from triangular matrix

                set(findobj('tag','hlasky'),'string',...
                        '','fontsize',8,'foregroundcolor','k'); drawnow;
  
                %Computation of commission error on a regular grid    
                if volbagridcheck==1      
                    
                    %Entered coordinates of the grid
                    if GUI==1 %If working with the GUI
                        fimin=str2num(get(findobj('tag','fimin'),'string')); 
                        fistep=str2num(get(findobj('tag','fistep'),'string'));
                        fimax=str2num(get(findobj('tag','fimax'),'string'));
                        lambdamin=str2num(get(findobj('tag','lambdamin'),'string')); 
                        lambdastep=str2num(get(findobj('tag','lambdastep'),'string'));
                        lambdamax=str2num(get(findobj('tag','lambdamax'),'string'));
                        h=str2num(get(findobj('tag','hgrid'),'string'));
                    elseif GUI==0 %If working without the GUI
                        fimin=lat_min; 
                        fistep=lat_step;
                        fimax=lat_max;
                        lambdamin=lon_min; 
                        lambdastep=lon_step;
                        lambdamax=lon_max;
                        clear lat_min lat_step lat_max lon_min lon_step lon_max
                    end
                    
                    %Input coordinates and some error checks
                    if GUI==0 && strcmp(fistep,'empty') && strcmp(fimax,'empty')
                        
                        if isempty(fimin) || ~isnumeric(fimin) || ~isvector(fimin)
                            error('When using a vector to define latitudes, the variable "lat_min" must be a real-valued scalar of vector.')
                        end
                        
                        %Vector of latitudes defined via the variable
                        %"lat_min"
                        fi=fimin(:); 
                        
                    else
                        if isempty(fimin) || ~isnumeric(fimin) || length(fimin)>1
                            if GUI==1
                                errordlg('The value of "Lat. min" must be a real-valued scalar.',...
                                    'Error in irregular surface selection panel');
                                error('The value of "Lat. min" must be a real-valued scalar.');
                            elseif GUI==0
                                error('The variable "lat_min" must be a real-valued scalar.');
                            end
                        end
                        if isempty(fistep) || ~isnumeric(fistep) || length(fistep)>1
                            if GUI==1
                                errordlg('The value of "Lat. step" must be a real-valued scalar.',...
                                    'Error in irregular surface selection panel');
                                error('The value of "Lat. step" must be a real-valued scalar.');
                            elseif GUI==0
                                error('The variable "lat_step" must be a real-valued scalar.');
                            end
                        end
                        if isempty(fimax) || ~isnumeric(fimax) || length(fimax)>1
                            if GUI==1
                                errordlg('The value of "Lat. max" must be a real-valued scalar.',...
                                    'Error in irregular surface selection panel');
                                error('The value of "Lat. max" must be a real-valued scalar.');
                            elseif GUI==0
                                error('The variable "lat_max" must be a real-valued scalar.');
                            end
                        end

                        if fimin>fimax
                            if GUI==1
                                errordlg('The value of "Lat. min" must be smaller than the "Lat. max" value.',...
                                    'Error in irregular surface selection panel');
                                error('The value of "Lat. min" must be smaller than the "Lat. max" value.'); 
                            elseif GUI==0
                                error('The value of "lat_min" must be smaller than the "lat_max" value.'); 
                            end
                        end
                        if fistep<=0
                            if GUI==1
                                errordlg('The value of "Lat. step" must be larger than zero.',...
                                    'Error in irregular surface selection panel');
                                error('The value of "Lat. step" must be larger than zero.');
                            elseif GUI==0
                                error('The value of "lat_step" must be larger than zero.');
                            end
                        end

                        if fimin>90 || fimin<-90
                            if GUI==1
                                errordlg('The value of "Lat. min" must be within the interval <-90 deg, 90 deg>.',...
                                    'Error in irregular surface selection panel');
                                error('The value of "Lat. min" must be within the interval <-90 deg, 90 deg>.');
                            elseif GUI==0
                                error('The value of "lat_min" must be within the interval <-90 deg, 90 deg>.');
                            end
                        end
                        if fimax>90 || fimax<-90
                            if GUI==1
                                errordlg('The value of "Lat. max" must be within the interval <-90 deg, 90 deg>.',...
                                    'Error in irregular surface selection panel');
                                error('The value of "Lat. max" must be within the interval <-90 deg, 90 deg>.');
                            elseif GUI==0
                                error('The value of "lat_max" must be within the interval <-90 deg, 90 deg>.');
                            end
                        end

                        fi=(fimin:fistep:fimax)';
                    end              
                    
                    if GUI==0 && strcmp(lambdastep,'empty') && strcmp(lambdamax,'empty')
                        
                        if isempty(lambdamin) || ~isnumeric(lambdamin) || ~isvector(lambdamin)
                            error('When using a vector to define longitudes, the variable "lon_min" must be a real-valued scalar of vector.')
                        end
                        
                        %Vector of latitudes defined via the variable
                        %"lat_min"
                        lambda=lambdamin(:); 
                        
                    else
                        if isempty(lambdamin) || ~isnumeric(lambdamin) || length(lambdamin)>1
                            if GUI==1
                                errordlg('The value of "Lon. min" must be a real-valued scalar.',...
                                    'Error in irregular surface selection panel');
                                error('The value of "Lon. min" must be a real-valued scalar.');
                            elseif GUI==0
                                error('The variable "lon_min" must be a real-valued scalar.');
                            end
                        end
                        if isempty(lambdastep) || ~isnumeric(lambdastep) || length(lambdastep)>1
                            if GUI==1
                                errordlg('The value of "Lon. step" must be a real-valued scalar.',...
                                    'Error in irregular surface selection panel');
                                error('The value of "Lon. step" must be a real-valued scalar.');
                            elseif GUI==0
                                error('The variable "lon_step" must be a real-valued scalar.');
                            end
                        end
                        if isempty(lambdamax) || ~isnumeric(lambdamax) || length(lambdamax)>1
                            if GUI==1
                                errordlg('The value of "Lon. max" must be a real-valued scalar.',...
                                    'Error in irregular surface selection panel');
                                error('The value of "Lon. max" must be a real-valued scalar.');
                            elseif GUI==0
                                error('The variable "lon_max" must be a real-valued scalar.');
                            end
                        end

                        if lambdamin>lambdamax
                            if GUI==1
                                errordlg('The value of "Lon. min" must be smaller than Lon. max value.',...
                                    'Error in irregular surface selection panel');
                                error('The value of "Lon. min" must be smaller than Lon. max value.');
                            elseif GUI==0
                                error('The value of "lon_min" must be smaller than Lon. max value.');
                            end
                        end
                        if lambdastep<=0
                            if GUI==1
                                errordlg('The value of "Lon. step" must be larger than zero.',...
                                    'Error in irregular surface selection panel');
                                error('The value of "Lon. step" must be larger than zero.');
                            elseif GUI==0
                                error('The value of "lon_step" must be larger than zero.');
                            end
                        end

                        if lambdamin>360 || lambdamin<-180
                            if GUI==1
                                errordlg('The value of "Lon. min" must be within the interval <-180 deg, 180 deg> or <0 deg, 360 deg>.',...
                                    'Error in irregular surface selection panel');
                                error('The value of "Lon. min" must be within the interval <-180 deg, 180 deg> or <0 deg, 360 deg>.');
                            elseif GUI==0
                                error('The value of "lon_min" must be within the interval <-180 deg, 180 deg> or <0 deg, 360 deg>.');
                            end
                        end
                        if lambdamax>360 || lambdamax<-180
                            if GUI==1
                                errordlg('The value of "Lon. max" must be within the interval <-180 deg, 180 deg> or <0 deg, 360 deg>.',...
                                    'Error in irregular surface selection panel');
                                error('The value of "Lon. max" must be within the interval <-180 deg, 180 deg> or <0 deg, 360 deg>.');
                            elseif GUI==0
                                error('The value of "lon_max" must be within the interval <-180 deg, 180 deg> or <0 deg, 360 deg>.');
                            end
                        end
                        if (lambdamax-lambdamin)>360
                            if GUI==1
                                errordlg('The longitudes must be in the range <-180 deg, 180 deg> or <0 deg, 360 deg>.',...
                                    'Error in irregular surface selection panel');
                                error('The longitudes must be in the range <-180 deg, 180 deg> or <0 deg, 360 deg>.');
                            elseif GUI==0
                                error('The longitudes must be in the range <-180 deg, 180 deg> or <0 deg, 360 deg>.');
                            end
                        end
                        
                        lambda=(lambdamin:lambdastep:lambdamax)'; 
                    end 
                    
                    if isempty(h) || ~isnumeric(h) || length(h)>1
                        if GUI==1
                            errordlg('The value of "Height above the reference surface" must be a real-valued scalar.',...
                                'Error in irregular surface selection panel');
                            error('The value of "Height above the reference surface" must be a real-valued scalar.');
                        elseif GUI==0
                            error('The variable "h" must be a real-valued scalar.');
                        end
                    end
                    
                    length_fi_plot=length(fi);
                    length_lambda=length(lambda);
                    
                    [fi,lambda]=meshgrid(fi,lambda);
                    [i_fi,j_fi]=size(fi);
                    [i_lambda,j_lambda]=size(lambda);
                    fi=fi'; 
                    fi=pi/180*(fi(:));
                    lambda=lambda'; 
                    lambda=pi/180*(lambda(:));  
                    
                    if coord==1
                        %Spherical radius
                        r=(R+h)*ones(length(fi),1);
                        hsph=h;
                        
                        %Spherical latitude
                        fiG=fi;  
                        
                        %Transform spherical latitude into the ellipsoidal
                        %latitude
                        [X,Y,Z]=sph2cart(0*fiG,fiG,r);
                        [fi, lambda_del, h]=cart2ell(X,Y,Z,[aEl eEl]);
                        
						clear X Y Z lambda_del
                    elseif coord==0
                        %Transformation of the ellipsoidal coordinates into the
                        %cartesian coordinates
                        [X,Y,Z]=ell2cart(fi,lambda,h,[aEl eEl]);
                        r=sqrt(X.*X+Y.*Y+Z.*Z); %Spherical radius
                        fiG=atan(Z./sqrt(X.*X+Y.*Y)); %Spherical latitude

                        clear X Y Z
                    end
                    
                    %If geoid/quasigeoid is to be computed
                    if any(volbapar==10) || any(volbapar==23)
                        if any(h~=0)
                            if GUI==1 %If working with the GUI
                                errordlg('In order to compute commission error of Geoid_undulation or Height_anomaly, the ellipsoidal height must be set to zero.',...
                                    'Error in point type selection');
                            end
                            error('In order to compute commission error of Geoid_undulation or Height_anomaly, the ellipsoidal height must be set to zero.');
                        end
                    end
                end
                
                %Computation of commission error in point-wise                            
                if volbadiskcheck==1               
                    %Ellipsoidal coordinates
                    if GUI==1 %If working with the GUI
                        fi=str2num(get(findobj('tag','fi'),'string'))';
                        lambda=str2num(get(findobj('tag','lambda'),'string'))';
                        h=str2num(get(findobj('tag','hdisk'),'string'))';
                    elseif GUI==0 %If working without the GUI
                        fi=lat(:);
                        lambda=lon(:);
                        h=h2(:);
                        clear lat lon h2
                    end
                    
                    length_fi=length(fi);
                    length_lambda=length(lambda);
                end

                %Computation of commission error from imported data 
                if volbaloadcheck==1
                    
                    if ~exist([loadadresar,loadname],'file') %Check whether the input data point file exists
                        if GUI==1 %Working with the GUI
                            errordlg('The entered data point file does not exist.',...
                                'Error in point type selection');
                        end
                        error('The entered data point file does not exist.')
                    end
                    if strcmp(loadname(end-3:end),'.mat') %Loading MAT file
                        Import=load([loadadresar,loadname]);
                        Import=struct2cell(Import);
                        Import=cell2mat(Import);
                        Import=Import(:,1:3);
                    else
                        Import=load([loadadresar,loadname]);
                    end
                    
                    %Ellipsoidal coordinates
                    fi=Import(:,1);
                    lambda=Import(:,2);
                    h=Import(:,3);                    
                    
                    volbadiskcheck=1;
                    
                    length_fi=length(fi);
                    length_lambda=length(lambda);
                end

                if volbaloadcheck==1 || volbadiskcheck==1
                    
                    if any(fi>90) || any(fi<-90)
                        if GUI==1 %If working with the GUI
                            errordlg('The values of latitude must be within the interval <-90 deg,90 deg>.',...
                                'Error in point type selection');
                            error('The values of latitude must be within the interval <-90 deg, 90 deg>.');
                        elseif GUI==0
                            error('The values of the "lat" variable must be within the interval <-90 deg, 90 deg>.');
                        end
                    end
                    if any(lambda>360) || any(lambda<-180)
                        if GUI==1 %If working with the GUI
                            errordlg('The values of longitudes must be within the interval <-180 deg, 180 deg> or <0 deg, 360 deg>.',...
                                'Error in point type selection');
                            error('The values of longitudes must be within the interval <-180 deg, 180 deg> or <0 deg, 360 deg>.');
                        elseif GUI==0
                            error('The values of the "lon" variable must be within the interval <-180 deg, 180 deg> or <0 deg, 360 deg>.');
                        end
                    end
                    
                    fi=pi/180*(fi(:));
                    lambda=pi/180*(lambda(:));  
                    
                    if coord==1 %Entered spherical coordinates                       
                        %Spherical radius
                        r=h;

                        %Spherical latitude
                        fiG=fi;  
                        
                        %Transform spherical latitude into the ellipsoidal
                        %latitude
                        [X,Y,Z]=sph2cart(lambda,fiG,r);
                        [fi,lambda_del,h]=cart2ell(X,Y,Z,[aEl eEl]);

                        clear X Y Z lambda_del
                    elseif coord==0 %Entered ellipsoidal coordinates
                        %Trasformation of (fi, lambda, h) into (X, Y, Z)
                        [X,Y,Z]=ell2cart(fi,lambda,h,[aEl eEl]);  
                        r=sqrt(X.*X+Y.*Y+Z.*Z); %Radius

                        %Spherical latitude
                        fiG=atan(Z./sqrt(X.*X+Y.*Y)); 
                        
                        clear X Y Z 
                    end
                end
                length_fiG=length(fiG);
                %=============================================================
                if volbadiskcheck==1                             

                    %Error message for input ellipsoidal coordinates
                    if isempty(fi)
                        if GUI==1
                            errordlg('The (ellipsoidal/spherical) "Latitude" array cannot be empty.',...
                                'Error in point type selection');
                            error('The (ellipsoidal/spherical) "Latitude" array cannot be empty.');
                        elseif GUI==0
                            error('The "lat" variable cannot be empty.');
                        end
                    end
                    if isempty(lambda)
                        if GUI==1
                            errordlg('The "Longitude" array cannot be empty.',...
                                'Error in point type selection');
                            error('The "Longitude" array cannot be empty.');
                        elseif GUI==0
                            error('The "lon" variable cannot be empty.');
                        end
                    end
                    if isempty(h)
                        if GUI==1
                            errordlg('The "Ellipsoidal height/Spherical radius" array cannot be empty.',...
                                'Error in point type selection');
                            error('The "Ellipsoidal height/Spherical radius" array cannot be empty.');
                        elseif GUI==0
                            error('The "h2" variable cannot be empty.');
                        end
                    end

                    if length_fi~=length_lambda || length_fi~=length(h) || length_lambda~=length(h)
                        if GUI==1 %If working with the GUI
                            errordlg('Ellipsoidal coordinates dimensions are not consistent or incorrectly entered.',...
                                'Error in point type selection')
                        end
                        error('Ellipsoidal coordinates dimensions are not consistent or incorrectly entered.')     
                    end               
                    
                    %If geoid/quasigeoid is to be computed
                    if any(volbapar==10) || any(volbapar==23)
                        if any(h~=0)
                            if GUI==1 %If working with the GUI
                                errordlg('In order to compute commission error of Geoid_undulation or Height_anomaly, the ellipsoidal height must be set to zero.',...
                                    'Error in point type selection');
                            end
                            error('In order to compute commission error of Geoid_undulation or Height_anomaly, the ellipsoidal height must be set to zero.');
                        end
                    end
                end
                
                %Initialization of the matrices and vectors for the computation of fnALFs
                length_fi=length(fi);
                length_fiG=length(fiG);
                Pnm=zeros(length_fi,nmax+1);
                q=(R./r);
                q2=(R./r).^2;
                u=cos(fiG);
                t=sin(fiG);
                
                %Initialization of the matrices and vectors for the 
                %computation of the first-order derivatives of fnALFs
                if any(volbapar==3) || any(volbapar==4) || any(volbapar==6) || any(volbapar==7) || any(volbapar==12) || any(volbapar==13) || any(volbapar==16)  || any(volbapar==20)
                    dALFs=1;
                    dPnm=zeros(length_fi,nmax+1);
                    qu=q./u;
                    tu=t./u;
                    
                    %Treatment of the dPnm singularity
                    singdPnm=fi==pi/2 | fi==-pi/2;
                else
                    dALFs=0;
                end 
                
                %Initialization of the matrices and vectors for the 
                %computation of the second-order derivatives of fnALFs
                if any(volbapar==6) || any(volbapar==12)
                    ddALFs=1;
                    ddPnm=zeros(length_fi,nmax+1);
                    
                    %Treatment of the ddPnm singularity
                    singddPnm=fi==pi/2 | fi==-pi/2;
                else
                    ddALFs=0;
                end  
                
                z=((nmax+1)*(nmax+2)-6)/2-(nmax-1);
                k=z;
                
                %Status line
                progressbar=findobj('tag','hlasky');
                
                if GUI==0 && Status_bar==1 %If working without the GUI
                    fprintf('Progress: m = ')
                end
                
                %% Summation over m
                for m=nmax:-1:0
                    
                    %Update of the progress bar
                    if GUI==1 && rem(m,10)==0 %If working with the GUI
                        set(progressbar,'string',...
                            sprintf('Progress: m = %5.0d',m),...
                            'fontsize',8); drawnow;
                    end
                    
                    if GUI==0 && Status_bar==1 %If working without the GUI
                        if rem(m,10)==0
                            fprintf('%d,',m)
                        end
                    end
                    
                    if m==1
                        z=z+1;
                    end
                    
                    %Standard forward column method
                    if m==0
                        Pnm(:,1)=1;
                    elseif m==1                    
                        Pnm(:,1)=sqrt(3)*u.*q;  
                    elseif m>1                            
                        i=2*(2:m);
                        i1=sqrt((i+ones(size(i)))./i);
                        Pnm(:,1)=u.^m*sqrt(3)*prod(i1).*q.^m;
                    end

                    if m==nmax
                    elseif m<=(nmax-1)
                        n=m+1;
                        anm=sqrt((2*n-1)*(2*n+1)/((n-m)*(n+m)));
                        Pnm(:,2)=anm*t.*Pnm(:,1).*q;
                    end

                    if m<(nmax-1)
                        j=3;
                        for n=m+2:nmax                            
                            anm=sqrt((2*n-1)*(2*n+1)/((n-m)*(n+m)));
                            bnm=sqrt((2*n+1)*(n+m-1)*(n-m-1)/((n-m)*(n+m)*(2*n-3)));
                            Pnm(:,j)=anm*t.*Pnm(:,j-1).*q-bnm*Pnm(:,j-2).*q2;
                            j=j+1;
                        end
                    end               
                      
                    %% Computation of the first-order derivatives of modified fnALFs
                    if dALFs==1  
                        if volbaALFs==1 || volbaALFs==3
                            enm=sqrt((2*(m:nmax)+1).*((m:nmax)-m).*((m:nmax)+m)./(2*(m:nmax)-1));

                            if m==0 %Zonal modified dALFs
                                dPnm(:,1)=0.*u;
                                dPnm(:,2)=sqrt(3)*u.*q;
                                dPnm(:,3:end)=-bsxfun(@times,2:nmax,tu).*Pnm(:,3:end)+(bsxfun(@times,enm(3:end),qu)).*Pnm(:,2:(end-1));
                            elseif m==nmax %Sectorial modified dALFs
                                dPnm(:,1)=-m*(tu).*Pnm(:,1);
                            else
                                dPnm(:,1)=-m*(tu).*Pnm(:,1); %Tesseral modified dALFs
                                dPnm(:,2:(nmax-m+1))=-bsxfun(@times,(m+1):nmax,tu).*Pnm(:,2:(nmax-m+1))+bsxfun(@times,enm(2:end),qu).*Pnm(:,1:(nmax-m)); %Teseralne ALFs
                            end

                        elseif volbaALFs==2
                            enm=sqrt((2*(m:nmax)+1).*((m:nmax)-m).*((m:nmax)+m)./(2*(m:nmax)-1));

                            if m==0 %Zonal modified dALFs
                                dPnm(:,1)=0.*u*1e-280;
                                dPnm(:,2)=sqrt(3)*u.*q*1e-280;
                                dPnm(:,3:end)=-bsxfun(@times,2:nmax,tu).*Pnm(:,3:end)+(bsxfun(@times,enm(3:end),qu)).*Pnm(:,2:(end-1));
                            elseif m==nmax
                                dPnm(:,1)=-m*(tu).*Pnm(:,1); %Sectorial modified dALFs
                            else
                                dPnm(:,1)=-m*(tu).*Pnm(:,1); %Tesseral modified dALFs
                                dPnm(:,2:(nmax-m+1))=-bsxfun(@times,(m+1):nmax,tu).*Pnm(:,2:(nmax-m+1))+bsxfun(@times,enm(2:end),qu).*Pnm(:,1:(nmax-m)); %Teseralne dALFs
                            end
                        end                 

                        %Treatment of the dPnm singularity
                        dPnm(singdPnm,:)=0;
                        
                        
                        %Computation of the second-order derivatives of modified fnALFs
                        if ddALFs==1
                            if m==0 %Zonal modified ddALFs
                                ddPnm=bsxfun(@times,tu,dPnm)-bsxfun(@times,(0:nmax).*((0:nmax)+1),Pnm);
                            else                                  
                                ddPnm(:,1:end-m)=bsxfun(@times,tu,dPnm(:,1:end-m))+bsxfun(@times,m^2./u.^2,Pnm(:,1:end-m))-bsxfun(@times,(m:nmax).*((m:nmax)+1),Pnm(:,1:end-m));
                            end
                                                                
                            %Treatment of the ddPnm singularity
                            ddPnm(singddPnm,:)=0;                    
                        end
                        
                    end
                    
                    %If nmin>2
                    if nmin~=2
                        Pnm(:,1:(nmin-m))=0;
                        
                        if dALFs==1
                            dPnm(:,1:(nmin-m))=0;
                            if ddALFs==1
                                ddPnm(:,1:(nmin-m))=0;
                            end
                        end
                    end
             
                    cosla=cos(m*lambda);
                    sinla=sin(m*lambda);

                    %% Loop for 1:NF (number of computing functionals)
                    for i=1:pocetpar                   
                        if volbapar(i)==1   
                        elseif volbapar(i)==2 %Deflection of the vertical eta                           

                            if m==nmax
                                ACeta=zeros(length_fiG,((nmax+1)*(nmax+2)-6)/2-(nmax-1));
                                ASeta=ACeta;
                                Aeta=zeros(length_fiG,cols_cov-(nmax-1));
                            end

                            if m==0
                                Aeta(:,1:2:end)=ACeta;
                                Aeta(:,2:2:end)=ASeta;

                                Aeta=[bsxfun(@times,m*Pnm(:,(3-m):(end-m)),sinla) Aeta];
                            elseif m==1
                                ACeta(:,z:k)=bsxfun(@times,m*Pnm(:,(3-m):(end-m)),sinla);
                                ASeta(:,z:k)=-bsxfun(@times,m*Pnm(:,(3-m):(end-m)),cosla); 
                            else
                                ACeta(:,z:k)=bsxfun(@times,m*Pnm(:,1:(nmax+1-m)),sinla);
                                ASeta(:,z:k)=-bsxfun(@times,m*Pnm(:,1:(nmax+1-m)),cosla);
                            end

                            if m==0
                                clear ACeta ASeta
                            end

                        elseif volbapar(i)==3 %Deflection of the vertical xi

                            if m==nmax
                                ACksi=zeros(length_fiG,((nmax+1)*(nmax+2)-6)/2-(nmax-1));
                                ASksi=ACksi;
                                Aksi=zeros(length_fiG,cols_cov-(nmax-1));
                            end

                            if m==0
                                Aksi(:,1:2:end)=ACksi;
                                Aksi(:,2:2:end)=ASksi;

                                Aksi=[bsxfun(@times,dPnm(:,(3-m):(end-m)),cosla) Aksi];
                            elseif m==1
                                ACksi(:,z:k)=bsxfun(@times,dPnm(:,(3-m):(end-m)),cosla);
                                ASksi(:,z:k)=bsxfun(@times,dPnm(:,(3-m):(end-m)),sinla); 
                            else
                                ACksi(:,z:k)=bsxfun(@times,dPnm(:,1:(nmax+1-m)),cosla);
                                ASksi(:,z:k)=bsxfun(@times,dPnm(:,1:(nmax+1-m)),sinla);
                            end

                            if m==0
                                clear ACksi ASksi
                            end

                        elseif volbapar(i)==4 %Deflection of the vertical Theta

                            if m==nmax
                                ACTeta=zeros(length_fiG,((nmax+1)*(nmax+2)-6)/2-(nmax-1));
                                ASTeta=ACTeta;
                                ACTksi=ACTeta;
                                ASTksi=ACTeta;
                                ATeta=zeros(length_fiG,cols_cov-(nmax-1));
                                ATksi=ATeta;
                            end

                            if m==0
                                ATeta(:,1:2:end)=ACTeta;
                                ATeta(:,2:2:end)=ASTeta;
                                ATksi(:,1:2:end)=ACTksi;
                                ATksi(:,2:2:end)=ASTksi;

                                ATksi=[bsxfun(@times,dPnm(:,(3-m):(end-m)),cosla) ATksi];
                                ATeta=[bsxfun(@times,m*Pnm(:,(3-m):(end-m)),sinla) ATeta];
                            elseif m==1
                                ACTksi(:,z:k)=bsxfun(@times,dPnm(:,(3-m):(end-m)),cosla);
                                ASTksi(:,z:k)=bsxfun(@times,dPnm(:,(3-m):(end-m)),sinla); 
                                
                                ACTeta(:,z:k)=bsxfun(@times,m*Pnm(:,(3-m):(end-m)),sinla);
                                ASTeta(:,z:k)=-bsxfun(@times,m*Pnm(:,(3-m):(end-m)),cosla); 
                            else
                                ACTksi(:,z:k)=bsxfun(@times,dPnm(:,1:(nmax+1-m)),cosla);
                                ASTksi(:,z:k)=bsxfun(@times,dPnm(:,1:(nmax+1-m)),sinla);
                                
                                ACTeta(:,z:k)=bsxfun(@times,m*Pnm(:,1:(nmax+1-m)),sinla);
                                ASTeta(:,z:k)=-bsxfun(@times,m*Pnm(:,1:(nmax+1-m)),cosla);
                            end

                            if m==0
                                clear ACTeta ASTeta ACTksi ASTksi
                            end

                        elseif volbapar(i)==5 %Disturbing potential

                            if m==nmax
                                ACT=zeros(length_fiG,((nmax+1)*(nmax+2)-6)/2-(nmax-1));
                                AST=ACT;
                                AT=zeros(length_fiG,cols_cov-(nmax-1));
                            end

                            if m==0
                                AT(:,1:2:end)=ACT;
                                AT(:,2:2:end)=AST;

                                AT=[bsxfun(@times,Pnm(:,(3-m):(end-m)),cosla) AT];
                            elseif m==1
                                ACT(:,z:k)=bsxfun(@times,Pnm(:,(3-m):(end-m)),cosla);
                                AST(:,z:k)=bsxfun(@times,Pnm(:,(3-m):(end-m)),sinla); 
                            else
                                ACT(:,z:k)=bsxfun(@times,Pnm(:,1:(nmax+1-m)),cosla);
                                AST(:,z:k)=bsxfun(@times,Pnm(:,1:(nmax+1-m)),sinla);
                            end

                            if m==0
                                clear ACT AST
                            end
                            
                        elseif volbapar(i)==6 %Disturbing tensor Trr_Tpp_Tll
                            
                            if m==nmax
                                ACTrr=zeros(length_fiG,((nmax+1)*(nmax+2)-6)/2-(nmax-1));
                                ASTrr=ACTrr;
                                ACTfifi=ACTrr;
                                ASTfifi=ACTrr;
                                ATrr=zeros(length_fiG,cols_cov-(nmax-1));
                                ATfifi=ATrr;
                                ACTll=ACTrr;
                                ASTll=ACTrr;
                                ATll=ATrr;
                                ampl_Trr=zeros(length_fiG,nmax-1); %Damping factor
                                for n=2:nmax
                                    ampl_Trr(:,n-1)=(n+1)*(n+2);
                                end
                            end

                            if m==0
                                ATrr(:,1:2:end)=ACTrr;
                                ATrr(:,2:2:end)=ASTrr;

                                ATrr=[bsxfun(@times,Pnm(:,(3-m):(end-m)).*ampl_Trr,cosla) ATrr];
                                
                                ATfifi(:,1:2:end)=ACTfifi;
                                ATfifi(:,2:2:end)=ASTfifi;
                                
                                ATfifi=[bsxfun(@times,ddPnm(:,(3-m):(end-m)),cosla) ATfifi];
                                
                                ATll(:,1:2:end)=ACTll;
                                ATll(:,2:2:end)=ASTll;

                                ATll=[bsxfun(@times,m^2*Pnm(:,(3-m):(end-m)),cosla) ATll];
                            elseif m==1
                                ACTrr(:,z:k)=bsxfun(@times,Pnm(:,(3-m):(end-m)).*ampl_Trr,cosla);
                                ASTrr(:,z:k)=bsxfun(@times,Pnm(:,(3-m):(end-m)).*ampl_Trr,sinla); 
                                
                                ACTfifi(:,z:k)=bsxfun(@times,ddPnm(:,(3-m):(end-m)),cosla);
                                ASTfifi(:,z:k)=bsxfun(@times,ddPnm(:,(3-m):(end-m)),sinla); 
                                
                                ACTll(:,z:k)=bsxfun(@times,m^2*Pnm(:,(3-m):(end-m)),cosla);
                                ASTll(:,z:k)=bsxfun(@times,m^2*Pnm(:,(3-m):(end-m)),sinla); 
                            else
                                ACTrr(:,z:k)=bsxfun(@times,Pnm(:,1:(nmax+1-m)).*ampl_Trr(:,(m-1):end),cosla);
                                ASTrr(:,z:k)=bsxfun(@times,Pnm(:,1:(nmax+1-m)).*ampl_Trr(:,(m-1):end),sinla);
                                
                                ACTfifi(:,z:k)=bsxfun(@times,ddPnm(:,1:(nmax+1-m)),cosla);
                                ASTfifi(:,z:k)=bsxfun(@times,ddPnm(:,1:(nmax+1-m)),sinla);
                                
                                ACTll(:,z:k)=bsxfun(@times,m^2*Pnm(:,1:(nmax+1-m)),cosla);
                                ASTll(:,z:k)=bsxfun(@times,m^2*Pnm(:,1:(nmax+1-m)),sinla);
                            end

                            if m==0
                                clear ampl_Trr ACTrr ASTrr ACTfifi ...
                                    ASTfifi ACTll ASTll
                            end
                            
                        elseif volbapar(i)==7 %Disturbing tensor Trp_Trl_Tpl
                            
                            if m==nmax
                                ACTrfi=zeros(length_fiG,((nmax+1)*(nmax+2)-6)/2-(nmax-1));
                                ASTrfi=ACTrfi;
                                ATrfi=zeros(length_fiG,cols_cov-(nmax-1));
                                ACTrl=ACTrfi;
                                ASTrl=ACTrfi;
                                ATrl=ATrfi;
                                ACTfil=ACTrfi;
                                ASTfil=ACTrfi;
                                ATfil=ATrfi;
                                ampl_Tr=zeros(length_fiG,nmax-1); %Damping factor
                                for n=2:nmax
                                    ampl_Tr(:,n-1)=(n+1);
                                end
                            end

                            if m==0
                                ATrfi(:,1:2:end)=ACTrfi;
                                ATrfi(:,2:2:end)=ASTrfi;

                                ATrfi=[bsxfun(@times,dPnm(:,(3-m):(end-m)).*ampl_Tr,cosla) ATrfi];
                                
                                ATrl(:,1:2:end)=ACTrl;
                                ATrl(:,2:2:end)=ASTrl;

                                ATrl=[bsxfun(@times,m*Pnm(:,(3-m):(end-m)).*ampl_Tr,sinla) ATrl];
                                
                                ATfil(:,1:2:end)=ACTfil;
                                ATfil(:,2:2:end)=ASTfil;

                                ATfil=[bsxfun(@times,m*dPnm(:,(3-m):(end-m)),sinla) ATfil];
                            elseif m==1
                                ACTrfi(:,z:k)=bsxfun(@times,dPnm(:,(3-m):(end-m)).*ampl_Tr,cosla);
                                ASTrfi(:,z:k)=bsxfun(@times,dPnm(:,(3-m):(end-m)).*ampl_Tr,sinla); 

                                ACTrl(:,z:k)=bsxfun(@times,m*Pnm(:,(3-m):(end-m)).*ampl_Tr,sinla);
                                ASTrl(:,z:k)=bsxfun(@times,m*Pnm(:,(3-m):(end-m)).*ampl_Tr,cosla); 
                                
                                ACTfil(:,z:k)=bsxfun(@times,m*dPnm(:,(3-m):(end-m)),sinla);
                                ASTfil(:,z:k)=bsxfun(@times,m*dPnm(:,(3-m):(end-m)),cosla); 
                            else
                                ACTrfi(:,z:k)=bsxfun(@times,dPnm(:,1:(nmax+1-m)).*ampl_Tr(:,(m-1):end),cosla);
                                ASTrfi(:,z:k)=bsxfun(@times,dPnm(:,1:(nmax+1-m)).*ampl_Tr(:,(m-1):end),sinla);

                                ACTrl(:,z:k)=bsxfun(@times,m*Pnm(:,1:(nmax+1-m)).*ampl_Tr(:,(m-1):end),sinla);
                                ASTrl(:,z:k)=bsxfun(@times,m*Pnm(:,1:(nmax+1-m)).*ampl_Tr(:,(m-1):end),cosla);
                                
                                ACTfil(:,z:k)=bsxfun(@times,m*dPnm(:,1:(nmax+1-m)),sinla);
                                ASTfil(:,z:k)=bsxfun(@times,m*dPnm(:,1:(nmax+1-m)),cosla);
                            end

                            if m==0
                                clear ampl_Trfi ACTrfi ASTrfi ACTrl ASTrl ... 
                                    ACTfil ASTfil ampl_Tr
                            end
                            
                        elseif volbapar(i)==8 %Disturbing tensor Txx_Tyy_Tzz
                        elseif volbapar(i)==9 %Disturbing tensor Txy_Txz_Tyz
                        elseif volbapar(i)==10 %Geoid undulation

                            if m==nmax
                                ACN=zeros(length_fiG,((nmax+1)*(nmax+2)-6)/2-(nmax-1));
                                ASN=ACN;
                                AN=zeros(length_fiG,cols_cov-(nmax-1));
                            end

                            if m==0
                                AN(:,1:2:end)=ACN;
                                AN(:,2:2:end)=ASN;

                                AN=[bsxfun(@times,Pnm(:,(3-m):(end-m)),cosla) AN];
                            elseif m==1
                                ACN(:,z:k)=bsxfun(@times,Pnm(:,(3-m):(end-m)),cosla);
                                ASN(:,z:k)=bsxfun(@times,Pnm(:,(3-m):(end-m)),sinla); 
                            else
                                ACN(:,z:k)=bsxfun(@times,Pnm(:,1:(nmax+1-m)),cosla);
                                ASN(:,z:k)=bsxfun(@times,Pnm(:,1:(nmax+1-m)),sinla);
                            end

                            if m==0
                                clear ACN ASN
                            end
                            
                        elseif volbapar(i)==11 %Gravitational potential

                            if m==nmax
                                ACV=zeros(length_fiG,((nmax+1)*(nmax+2)-6)/2-(nmax-1));
                                ASV=ACV;
                                AV=zeros(length_fiG,cols_cov-(nmax-1));
                            end

                            if m==0
                                AV(:,1:2:end)=ACV;
                                AV(:,2:2:end)=ASV;

                                AV=[bsxfun(@times,Pnm(:,(3-m):(end-m)),cosla) AV];
                            elseif m==1
                                ACV(:,z:k)=bsxfun(@times,Pnm(:,(3-m):(end-m)),cosla);
                                ASV(:,z:k)=bsxfun(@times,Pnm(:,(3-m):(end-m)),sinla); 
                            else
                                ACV(:,z:k)=bsxfun(@times,Pnm(:,1:(nmax+1-m)),cosla);
                                ASV(:,z:k)=bsxfun(@times,Pnm(:,1:(nmax+1-m)),sinla);
                            end

                            if m==0
                                clear ACV ASV
                            end

                        elseif volbapar(i)==12 %Gravitational tensor Vrr_Vpp_Vll
                            
                            if m==nmax
                                ACVrr=zeros(length_fiG,((nmax+1)*(nmax+2)-6)/2-(nmax-1));
                                ASVrr=ACVrr;
                                ACVfifi=ACVrr;
                                ASVfifi=ACVrr;
                                AVrr=zeros(length_fiG,cols_cov-(nmax-1));
                                AVfifi=AVrr;
                                ACVll=ACVrr;
                                ASVll=ACVrr;
                                AVll=AVrr;
                                ampl_Vrr=zeros(length_fiG,nmax-1); %Damping factor
                                for n=2:nmax
                                    ampl_Vrr(:,n-1)=(n+1)*(n+2);
                                end
                            end

                            if m==0
                                AVrr(:,1:2:end)=ACVrr;
                                AVrr(:,2:2:end)=ASVrr;

                                AVrr=[bsxfun(@times,Pnm(:,(3-m):(end-m)).*ampl_Vrr,cosla) AVrr];
                                
                                AVfifi(:,1:2:end)=ACVfifi;
                                AVfifi(:,2:2:end)=ASVfifi;
                                
                                AVfifi=[bsxfun(@times,ddPnm(:,(3-m):(end-m)),cosla) AVfifi];
                                
                                AVll(:,1:2:end)=ACVll;
                                AVll(:,2:2:end)=ASVll;

                                AVll=[bsxfun(@times,m^2*Pnm(:,(3-m):(end-m)),cosla) AVll];
                            elseif m==1
                                ACVrr(:,z:k)=bsxfun(@times,Pnm(:,(3-m):(end-m)).*ampl_Vrr,cosla);
                                ASVrr(:,z:k)=bsxfun(@times,Pnm(:,(3-m):(end-m)).*ampl_Vrr,sinla); 
                                
                                ACVfifi(:,z:k)=bsxfun(@times,ddPnm(:,(3-m):(end-m)),cosla);
                                ASVfifi(:,z:k)=bsxfun(@times,ddPnm(:,(3-m):(end-m)),sinla); 
                                
                                ACVll(:,z:k)=bsxfun(@times,m^2*Pnm(:,(3-m):(end-m)),cosla);
                                ASVll(:,z:k)=bsxfun(@times,m^2*Pnm(:,(3-m):(end-m)),sinla); 
                            else
                                ACVrr(:,z:k)=bsxfun(@times,Pnm(:,1:(nmax+1-m)).*ampl_Vrr(:,(m-1):end),cosla);
                                ASVrr(:,z:k)=bsxfun(@times,Pnm(:,1:(nmax+1-m)).*ampl_Vrr(:,(m-1):end),sinla);
                                
                                ACVfifi(:,z:k)=bsxfun(@times,ddPnm(:,1:(nmax+1-m)),cosla);
                                ASVfifi(:,z:k)=bsxfun(@times,ddPnm(:,1:(nmax+1-m)),sinla);
                                
                                ACVll(:,z:k)=bsxfun(@times,m^2*Pnm(:,1:(nmax+1-m)),cosla);
                                ASVll(:,z:k)=bsxfun(@times,m^2*Pnm(:,1:(nmax+1-m)),sinla);
                            end

                            if m==0
                                clear ampl_Vrr ACVrr ASVrr ACVfifi ...
                                    ASVfifi ACVll ASVll
                            end
                            
                        elseif volbapar(i)==13 %Gravitational tensor Vrp_Vrl_Vpl
                            
                            if m==nmax
                                ACVrfi=zeros(length_fiG,((nmax+1)*(nmax+2)-6)/2-(nmax-1));
                                ASVrfi=ACVrfi;
                                AVrfi=zeros(length_fiG,cols_cov-(nmax-1));
                                ACVrl=ACVrfi;
                                ASVrl=ACVrfi;
                                AVrl=AVrfi;
                                ACVfil=ACVrfi;
                                ASVfil=ACVrfi;
                                AVfil=AVrfi;
                                ampl_Vr=zeros(length_fiG,nmax-1); %Damping factor
                                for n=2:nmax
                                    ampl_Vr(:,n-1)=(n+1);
                                end
                            end

                            if m==0
                                AVrfi(:,1:2:end)=ACVrfi;
                                AVrfi(:,2:2:end)=ASVrfi;

                                AVrfi=[bsxfun(@times,dPnm(:,(3-m):(end-m)).*ampl_Vr,cosla) AVrfi];
                                
                                AVrl(:,1:2:end)=ACVrl;
                                AVrl(:,2:2:end)=ASVrl;

                                AVrl=[bsxfun(@times,m*Pnm(:,(3-m):(end-m)).*ampl_Vr,sinla) AVrl];
                                
                                AVfil(:,1:2:end)=ACVfil;
                                AVfil(:,2:2:end)=ASVfil;

                                AVfil=[bsxfun(@times,m*dPnm(:,(3-m):(end-m)),sinla) AVfil];
                            elseif m==1
                                ACVrfi(:,z:k)=bsxfun(@times,dPnm(:,(3-m):(end-m)).*ampl_Vr,cosla);
                                ASVrfi(:,z:k)=bsxfun(@times,dPnm(:,(3-m):(end-m)).*ampl_Vr,sinla); 

                                ACVrl(:,z:k)=bsxfun(@times,m*Pnm(:,(3-m):(end-m)).*ampl_Vr,sinla);
                                ASVrl(:,z:k)=bsxfun(@times,m*Pnm(:,(3-m):(end-m)).*ampl_Vr,cosla); 
                                
                                ACVfil(:,z:k)=bsxfun(@times,m*dPnm(:,(3-m):(end-m)),sinla);
                                ASVfil(:,z:k)=bsxfun(@times,m*dPnm(:,(3-m):(end-m)),cosla); 
                            else
                                ACVrfi(:,z:k)=bsxfun(@times,dPnm(:,1:(nmax+1-m)).*ampl_Vr(:,(m-1):end),cosla);
                                ASVrfi(:,z:k)=bsxfun(@times,dPnm(:,1:(nmax+1-m)).*ampl_Vr(:,(m-1):end),sinla);

                                ACVrl(:,z:k)=bsxfun(@times,m*Pnm(:,1:(nmax+1-m)).*ampl_Vr(:,(m-1):end),sinla);
                                ASVrl(:,z:k)=bsxfun(@times,m*Pnm(:,1:(nmax+1-m)).*ampl_Vr(:,(m-1):end),cosla);
                                
                                ACVfil(:,z:k)=bsxfun(@times,m*dPnm(:,1:(nmax+1-m)),sinla);
                                ASVfil(:,z:k)=bsxfun(@times,m*dPnm(:,1:(nmax+1-m)),cosla);
                            end

                            if m==0
                                clear ampl_Vrfi ACVrfi ASVrfi ACVrl ...
                                    ASVrl ACVfil ASVfil ampl_Vr
                            end
                            
                        elseif volbapar(i)==14 %Gravitational tensor Vxx_Vyy_Vzz
                        elseif volbapar(i)==15 %Vxy_Vxz_Vyz
                        elseif volbapar(i)==16 %Gravity vector gX_gY_gZ

                            if m==nmax
                                ACWr=zeros(length_fiG,((nmax+1)*(nmax+2)-6)/2-(nmax-1));
                                ASWr=ACWr;
                                ACWfi=ACWr;
                                ASWfi=ACWr;
                                ACWlambda=ACWr;
                                ASWlambda=ACWr;
                                AWr=zeros(length_fiG,cols_cov-(nmax-1));
                                AWfi=AWr;
                                AWlambda=AWr;
                                ampl_Wr=zeros(length_fiG,nmax-1); %Damping factor
                                for n=2:nmax
                                    ampl_Wr(:,n-1)=(n+1);
                                end
                            end

                            if m==0
                                AWr(:,1:2:end)=ACWr;
                                AWr(:,2:2:end)=ASWr;
                                AWfi(:,1:2:end)=ACWfi;
                                AWfi(:,2:2:end)=ASWfi;
                                AWlambda(:,1:2:end)=ACWlambda;
                                AWlambda(:,2:2:end)=ASWlambda;

                                AWlambda=[bsxfun(@times,m*Pnm(:,(3-m):(end-m)),sinla) AWlambda];
                                AWfi=[bsxfun(@times,dPnm(:,(3-m):(end-m)),cosla) AWfi];
                                AWr=[bsxfun(@times,Pnm(:,(3-m):(end-m)).*ampl_Wr,cosla) AWr];
                            elseif m==1
                                ACWlambda(:,z:k)=bsxfun(@times,m*Pnm(:,(3-m):(end-m)),sinla);
                                ASWlambda(:,z:k)=-bsxfun(@times,m*Pnm(:,(3-m):(end-m)),cosla); 
                                
                                ACWfi(:,z:k)=bsxfun(@times,dPnm(:,(3-m):(end-m)),cosla);
                                ASWfi(:,z:k)=bsxfun(@times,dPnm(:,(3-m):(end-m)),sinla);
                                
                                ACWr(:,z:k)=bsxfun(@times,Pnm(:,(3-m):(end-m)).*ampl_Wr,cosla);
                                ASWr(:,z:k)=bsxfun(@times,Pnm(:,(3-m):(end-m)).*ampl_Wr,sinla); 
                            else
                                ACWlambda(:,z:k)=bsxfun(@times,m*Pnm(:,1:(nmax+1-m)),sinla);
                                ASWlambda(:,z:k)=-bsxfun(@times,m*Pnm(:,1:(nmax+1-m)),cosla);
                                
                                ACWfi(:,z:k)=bsxfun(@times,dPnm(:,1:(nmax+1-m)),cosla);
                                ASWfi(:,z:k)=bsxfun(@times,dPnm(:,1:(nmax+1-m)),sinla);
                                
                                ACWr(:,z:k)=bsxfun(@times,Pnm(:,1:(nmax+1-m)).*ampl_Wr(:,(m-1):end),cosla);
                                ASWr(:,z:k)=bsxfun(@times,Pnm(:,1:(nmax+1-m)).*ampl_Wr(:,(m-1):end),sinla);
                            end

                            if m==0
                                clear ampl_Wr ACWr ASWr ACWfi ASWfi ...
                                    ACWlambda ASWlambda
                            end
                                                        
                        elseif volbapar(i)==17 %Gravity sa

                            if m==nmax
                                ACg_sa=zeros(length_fiG,((nmax+1)*(nmax+2)-6)/2-(nmax-1));
                                ASg_sa=ACg_sa;
                                Ag_sa=zeros(length_fiG,cols_cov-(nmax-1));
                                ampl_g_sa=zeros(length_fiG,nmax-1); %Damping factor
                                for n=2:nmax
                                    ampl_g_sa(:,n-1)=(n+1);
                                end
                            end

                            if m==0
                                Ag_sa(:,1:2:end)=ACg_sa;
                                Ag_sa(:,2:2:end)=ASg_sa;

                                Ag_sa=[bsxfun(@times,Pnm(:,(3-m):(end-m)).*ampl_g_sa,cosla) Ag_sa];
                            elseif m==1
                                ACg_sa(:,z:k)=bsxfun(@times,Pnm(:,(3-m):(end-m)).*ampl_g_sa,cosla);
                                ASg_sa(:,z:k)=bsxfun(@times,Pnm(:,(3-m):(end-m)).*ampl_g_sa,sinla); 
                            else
                                ACg_sa(:,z:k)=bsxfun(@times,Pnm(:,1:(nmax+1-m)).*ampl_g_sa(:,(m-1):end),cosla);
                                ASg_sa(:,z:k)=bsxfun(@times,Pnm(:,1:(nmax+1-m)).*ampl_g_sa(:,(m-1):end),sinla);
                            end

                            if m==0
                                clear ampl_g_sa ACg_sa ASg_sa
                            end

                        elseif volbapar(i)==18 %Gravity potential                            

                            if m==nmax
                                ACW=zeros(length_fiG,((nmax+1)*(nmax+2)-6)/2-(nmax-1));
                                ASW=ACW;
                                AW=zeros(length_fiG,cols_cov-(nmax-1));
                            end

                            if m==0
                                AW(:,1:2:end)=ACW;
                                AW(:,2:2:end)=ASW;

                                AW=[bsxfun(@times,Pnm(:,(3-m):(end-m)),cosla) AW];
                            elseif m==1
                                ACW(:,z:k)=bsxfun(@times,Pnm(:,(3-m):(end-m)),cosla);
                                ASW(:,z:k)=bsxfun(@times,Pnm(:,(3-m):(end-m)),sinla); 
                            else
                                ACW(:,z:k)=bsxfun(@times,Pnm(:,1:(nmax+1-m)),cosla);
                                ASW(:,z:k)=bsxfun(@times,Pnm(:,1:(nmax+1-m)),sinla);
                            end

                            if m==0
                                clear ACW ASW
                            end

                        elseif volbapar(i)==19 %Gravity anomaly sa

                            if m==nmax
                                ACanomalia_sa=zeros(length_fiG,((nmax+1)*(nmax+2)-6)/2-(nmax-1));
                                ASanomalia_sa=ACanomalia_sa;
                                Aanomalia_sa=zeros(length_fiG,cols_cov-(nmax-1));
                                ampl_anomalia_sa=zeros(length_fiG,nmax-1); %Damping factor
                                for n=2:nmax
                                    ampl_anomalia_sa(:,n-1)=(n-1);
                                end
                            end

                            if m==0
                                Aanomalia_sa(:,1:2:end)=ACanomalia_sa;
                                Aanomalia_sa(:,2:2:end)=ASanomalia_sa;

                                Aanomalia_sa=[bsxfun(@times,Pnm(:,(3-m):(end-m)).*ampl_anomalia_sa,cosla) Aanomalia_sa];
                            elseif m==1
                                ACanomalia_sa(:,z:k)=bsxfun(@times,Pnm(:,(3-m):(end-m)).*ampl_anomalia_sa,cosla);
                                ASanomalia_sa(:,z:k)=bsxfun(@times,Pnm(:,(3-m):(end-m)).*ampl_anomalia_sa,sinla); 
                            else
                                ACanomalia_sa(:,z:k)=bsxfun(@times,Pnm(:,1:(nmax+1-m)).*ampl_anomalia_sa(:,(m-1):end),cosla);
                                ASanomalia_sa(:,z:k)=bsxfun(@times,Pnm(:,1:(nmax+1-m)).*ampl_anomalia_sa(:,(m-1):end),sinla);
                            end

                            if m==0
                                clear ampl_anomalia_sa ACanomalia_sa ...
                                    ASanomalia_sa
                            end

                        elseif volbapar(i)==20 %Gravity disturbance

                            if m==nmax
                                ACWr_porucha=zeros(length_fiG,((nmax+1)*(nmax+2)-6)/2-(nmax-1));
                                ASWr_porucha=ACWr_porucha;
                                ACWfi_porucha=ACWr_porucha;
                                ASWfi_porucha=ACWr_porucha;
                                ACWlambda_porucha=ACWr_porucha;
                                ASWlambda_porucha=ACWr_porucha;
                                AWr_porucha=zeros(length_fiG,cols_cov-(nmax-1));
                                AWfi_porucha=AWr_porucha;
                                AWlambda_porucha=AWr_porucha;
                                ampl_Wr_porucha=zeros(length_fiG,nmax-1); %Damping factor
                                for n=2:nmax
                                    ampl_Wr_porucha(:,n-1)=(n+1);
                                end
                            end

                            if m==0
                                AWr_porucha(:,1:2:end)=ACWr_porucha;
                                AWr_porucha(:,2:2:end)=ASWr_porucha;
                                AWfi_porucha(:,1:2:end)=ACWfi_porucha;
                                AWfi_porucha(:,2:2:end)=ASWfi_porucha;
                                AWlambda_porucha(:,1:2:end)=ACWlambda_porucha;
                                AWlambda_porucha(:,2:2:end)=ASWlambda_porucha;

                                AWlambda_porucha=[bsxfun(@times,m*Pnm(:,(3-m):(end-m)),sinla) AWlambda_porucha];
                                AWfi_porucha=[bsxfun(@times,dPnm(:,(3-m):(end-m)),cosla) AWfi_porucha];
                                AWr_porucha=[bsxfun(@times,Pnm(:,(3-m):(end-m)).*ampl_Wr_porucha,cosla) AWr_porucha];
                            elseif m==1
                                ACWlambda_porucha(:,z:k)=bsxfun(@times,m*Pnm(:,(3-m):(end-m)),sinla);
                                ASWlambda_porucha(:,z:k)=-bsxfun(@times,m*Pnm(:,(3-m):(end-m)),cosla); 
                                
                                ACWfi_porucha(:,z:k)=bsxfun(@times,dPnm(:,(3-m):(end-m)),cosla);
                                ASWfi_porucha(:,z:k)=bsxfun(@times,dPnm(:,(3-m):(end-m)),sinla);
                                
                                ACWr_porucha(:,z:k)=bsxfun(@times,Pnm(:,(3-m):(end-m)).*ampl_Wr_porucha,cosla);
                                ASWr_porucha(:,z:k)=bsxfun(@times,Pnm(:,(3-m):(end-m)).*ampl_Wr_porucha,sinla); 
                            else
                                ACWlambda_porucha(:,z:k)=bsxfun(@times,m*Pnm(:,1:(nmax+1-m)),sinla);
                                ASWlambda_porucha(:,z:k)=-bsxfun(@times,m*Pnm(:,1:(nmax+1-m)),cosla);
                                
                                ACWfi_porucha(:,z:k)=bsxfun(@times,dPnm(:,1:(nmax+1-m)),cosla);
                                ASWfi_porucha(:,z:k)=bsxfun(@times,dPnm(:,1:(nmax+1-m)),sinla);
                                
                                ACWr_porucha(:,z:k)=bsxfun(@times,Pnm(:,1:(nmax+1-m)).*ampl_Wr_porucha(:,(m-1):end),cosla);
                                ASWr_porucha(:,z:k)=bsxfun(@times,Pnm(:,1:(nmax+1-m)).*ampl_Wr_porucha(:,(m-1):end),sinla);
                            end

                            if m==0
                                clear ampl_Wr_porucha ACWr_porucha ...
                                    ASWr_porucha ACWfi_porucha ...
                                    ASWfi_porucha ACWlambda_porucha ASWlambda
                            end

                        elseif volbapar(i)==21 %Gravity disturbance sa

                            if m==nmax
                                ACporucha_sa=zeros(length_fiG,((nmax+1)*(nmax+2)-6)/2-(nmax-1));
                                ASporucha_sa=ACporucha_sa;
                                Aporucha_sa=zeros(length_fiG,cols_cov-(nmax-1));
                                ampl_porucha_sa=zeros(length_fiG,nmax-1); %Damping factor
                                for n=2:nmax
                                    ampl_porucha_sa(:,n-1)=(n+1);
                                end
                            end

                            if m==0
                                Aporucha_sa(:,1:2:end)=ACporucha_sa;
                                Aporucha_sa(:,2:2:end)=ASporucha_sa;

                                Aporucha_sa=[bsxfun(@times,Pnm(:,(3-m):(end-m)).*ampl_porucha_sa,cosla) Aporucha_sa];
                            elseif m==1
                                ACporucha_sa(:,z:k)=bsxfun(@times,Pnm(:,(3-m):(end-m)).*ampl_porucha_sa,cosla);
                                ASporucha_sa(:,z:k)=bsxfun(@times,Pnm(:,(3-m):(end-m)).*ampl_porucha_sa,sinla); 
                            else
                                ACporucha_sa(:,z:k)=bsxfun(@times,Pnm(:,1:(nmax+1-m)).*ampl_porucha_sa(:,(m-1):end),cosla);
                                ASporucha_sa(:,z:k)=bsxfun(@times,Pnm(:,1:(nmax+1-m)).*ampl_porucha_sa(:,(m-1):end),sinla);
                            end

                            if m==0
                                clear ampl_porucha_sa ACporucha_sa ASporucha_sa
                            end

                        elseif volbapar(i)==22 %Height anomaly Ell

                            if m==nmax
                                ACzetaEl=zeros(length_fiG,((nmax+1)*(nmax+2)-6)/2-(nmax-1));
                                ASzetaEl=ACzetaEl;
                                AzetaEl=zeros(length_fiG,cols_cov-(nmax-1));
                            end

                            if m==0
                                AzetaEl(:,1:2:end)=ACzetaEl;
                                AzetaEl(:,2:2:end)=ASzetaEl;

                                AzetaEl=[bsxfun(@times,Pnm(:,(3-m):(end-m)),cosla) AzetaEl];
                            elseif m==1
                                ACzetaEl(:,z:k)=bsxfun(@times,Pnm(:,(3-m):(end-m)),cosla);
                                ASzetaEl(:,z:k)=bsxfun(@times,Pnm(:,(3-m):(end-m)),sinla); 
                            else
                                ACzetaEl(:,z:k)=bsxfun(@times,Pnm(:,1:(nmax+1-m)),cosla);
                                ASzetaEl(:,z:k)=bsxfun(@times,Pnm(:,1:(nmax+1-m)),sinla);
                            end

                            if m==0
                                clear ACzetaEl ASzetaEl
                            end

                        elseif volbapar(i)==23 %Height anomaly
                            
                            if m==nmax
                                ACzeta=zeros(length_fiG,((nmax+1)*(nmax+2)-6)/2-(nmax-1));
                                ASzeta=ACzeta;
                                Azeta=zeros(length_fiG,cols_cov-(nmax-1));
                            end

                            if m==0
                                Azeta(:,1:2:end)=ACzeta;
                                Azeta(:,2:2:end)=ASzeta;

                                Azeta=[bsxfun(@times,Pnm(:,(3-m):(end-m)),cosla) Azeta];
                            elseif m==1
                                ACzeta(:,z:k)=bsxfun(@times,Pnm(:,(3-m):(end-m)),cosla);
                                ASzeta(:,z:k)=bsxfun(@times,Pnm(:,(3-m):(end-m)),sinla); 
                            else
                                ACzeta(:,z:k)=bsxfun(@times,Pnm(:,1:(nmax+1-m)),cosla);
                                ASzeta(:,z:k)=bsxfun(@times,Pnm(:,1:(nmax+1-m)),sinla);
                            end

                            if m==0
                                clear ACzeta ASzeta
                            end
                            
                        elseif volbapar(i)==24 %Second radial derivative of disturbing potential

                            if m==nmax
                                ACT_rr=zeros(length_fiG,((nmax+1)*(nmax+2)-6)/2-(nmax-1));
                                AST_rr=ACT_rr;
                                AT_rr=zeros(length_fiG,cols_cov-(nmax-1));
                                ampl_T_rr=zeros(length_fiG,nmax-1); %Damping factor
                                for n=2:nmax
                                    ampl_T_rr(:,n-1)=(n+1)*(n+2);
                                end
                            end

                            if m==0
                                AT_rr(:,1:2:end)=ACT_rr;
                                AT_rr(:,2:2:end)=AST_rr;

                                AT_rr=[bsxfun(@times,Pnm(:,(3-m):(end-m)).*ampl_T_rr,cosla) AT_rr];
                            elseif m==1
                                ACT_rr(:,z:k)=bsxfun(@times,Pnm(:,(3-m):(end-m)).*ampl_T_rr,cosla);
                                AST_rr(:,z:k)=bsxfun(@times,Pnm(:,(3-m):(end-m)).*ampl_T_rr,sinla); 
                            else
                                ACT_rr(:,z:k)=bsxfun(@times,Pnm(:,1:(nmax+1-m)).*ampl_T_rr(:,(m-1):end),cosla);
                                AST_rr(:,z:k)=bsxfun(@times,Pnm(:,1:(nmax+1-m)).*ampl_T_rr(:,(m-1):end),sinla);
                            end

                            if m==0
                                clear ampl_T_rr ACT_rr AST_rr
                            end

                        elseif volbapar(i)==25 %Second radial derivative of gravity potential

                            if m==nmax
                                ACWrr=zeros(length_fiG,((nmax+1)*(nmax+2)-6)/2-(nmax-1));
                                ASWrr=ACWrr;
                                AWrr=zeros(length_fiG,cols_cov-(nmax-1));
                                ampl_Wrr=zeros(length_fiG,nmax-1); %Damping factor
                                for n=2:nmax
                                    ampl_Wrr(:,n-1)=(n+1)*(n+2);
                                end
                            end

                            if m==0
                                AWrr(:,1:2:end)=ACWrr;
                                AWrr(:,2:2:end)=ASWrr;

                                AWrr=[bsxfun(@times,Pnm(:,(3-m):(end-m)).*ampl_Wrr,cosla) AWrr];
                            elseif m==1
                                ACWrr(:,z:k)=bsxfun(@times,Pnm(:,(3-m):(end-m)).*ampl_Wrr,cosla);
                                ASWrr(:,z:k)=bsxfun(@times,Pnm(:,(3-m):(end-m)).*ampl_Wrr,sinla); 
                            else
                                ACWrr(:,z:k)=bsxfun(@times,Pnm(:,1:(nmax+1-m)).*ampl_Wrr(:,(m-1):end),cosla);
                                ASWrr(:,z:k)=bsxfun(@times,Pnm(:,1:(nmax+1-m)).*ampl_Wrr(:,(m-1):end),sinla);
                            end

                            if m==0
                                clear ampl_Wrr ACWrr ASWrr
                            end
                        end
                    end
                    
                    if m==0
                    elseif m==1
                        z=z-nmax+m-1;
                        k=z+nmax-m; 
                    else
                        z=z-nmax+m-2;
                        k=z+nmax-m+1; 
                    end
                end
                
                %Update of progress bar
                if GUI==1 %If working with the GUI
                    set(progressbar,'string','Progress: Matrix multiplications...','fontsize',8); drawnow;
                end
                
                if GUI==0 && Status_bar==1 %If working without the GUI
                    fprintf('\n')
                    fprintf('Progress: Matrix multilications...\n')
                end
                
                clear Pnm dPnm ddPnm t u q q2 tu qu sinla cosla enm ...
                    singdPnm singddPnm
                
                %Computation of the normal gravity for eta, xi, Theta,
                %geoid undulation, zeta El, zeta
                if any(volbapar==2) || any(volbapar==3) || any(volbapar==4) || any(volbapar==10) || any(volbapar==22) || any(volbapar==23)                        
                    bEl=aEl*sqrt(1-eEl^2);
                    EEl=sqrt(aEl^2-bEl^2);
                    
                    if volbagridcheck==1
                        %Computation of the ellipsoidal harmonic coordinates
                        [X,Y,Z]=ell2cart(fi,0*zeros(length_fi,1),h,[aEl eEl]);
                        ugama=sqrt((X.*X+Y.*Y+Z.*Z-EEl.^2).*(1./2+1./2.*sqrt(1+((4.*EEl.^2.*Z.*Z)./(X.*X+Y.*Y+Z.*Z-EEl.^2).^2))));                        
                        betagama=atan((Z.*sqrt(ugama.^2+EEl.^2))./(ugama.*sqrt(X.*X+Y.*Y)));
                        clear X Y Z
                    elseif volbadiskcheck==1
                        %Computation of the ellipsoidal harmonic coordinates
                        [X,Y,Z]=ell2cart(fi,lambda,h,[aEl eEl]);
                        ugama=sqrt((X.*X+Y.*Y+Z.*Z-EEl.^2).*(1./2+1./2.*sqrt(1+((4.*EEl.^2.*Z.*Z)./(X.*X+Y.*Y+Z.*Z-EEl.^2).^2))));                        
                        betagama=atan((Z.*sqrt(ugama.^2+EEl.^2))./(ugama.*sqrt(X.*X+Y.*Y)));
                        clear X Y Z
                    end

                    wgama=sqrt((ugama.^2+EEl^2*sin(betagama).^2)./(ugama.^2+EEl^2));
                    qgama=1/2*((1+(3*ugama.^2)./EEl^2).*atan(EEl./ugama)-3*ugama./EEl);
                    qgama_=3*(1+(ugama.^2)./EEl^2).*(1-ugama./EEl.*atan(EEl./ugama))-1;
                    qgama0=1/2*((1+(3*bEl^2)/EEl^2)*atan(EEl/bEl)-3*bEl/EEl);
                        
                    gamau=-1./wgama.*(GMEl./(ugama.^2+EEl.^2)+(omegaEl.^2.*aEl.^2.*EEl)./(ugama.^2+EEl.^2).*(qgama_./qgama0).*(1./2.*sin(betagama).^2-1./6)-omegaEl.^2.*ugama.*cos(betagama).^2);
                    gamabeta=-1./wgama.*(-(omegaEl.^2.*aEl.^2.*qgama)./(sqrt(ugama.^2+EEl.^2).*qgama0)+omegaEl.^2.*sqrt(ugama.^2+EEl.^2)).*sin(betagama).*cos(betagama);
                        
                    clear ugama betagama wgama qgama qgama_ qgama0
                       
                    gamaP=sqrt(gamau.^2+gamabeta.^2);
                        
                    clear gamau gamabeta
                end
                
                %% Final computation of commission errors of the functionals             
                for i=1:pocetpar
                    if volbapar(i)==1                
                    elseif volbapar(i)==2 %Deflection of the vertical eta

                        Pg=sum((Aeta*covmat).*Aeta,2);
                        clear Aeta
                        Pg=(GM./(r.^2.*gamaP.*cos(fiG)).*sqrt(Pg))*(180/pi)*3600;
                        Pg(fi==pi/2 | fi==-pi/2)=0;
                        
                    elseif volbapar(i)==3 %Deflection of the vertical xi

                        Pg=sum((Aksi*covmat).*Aksi,2);
                        clear Aksi
                        Pg=(GM./(r.^2.*gamaP).*sqrt(Pg))*(180/pi)*3600;
                        
                    elseif volbapar(i)==4 %Deflection of the vertical Theta

                        eta=sum((ATeta*covmat).*ATeta,2);
                        clear ATeta
                       
                        ksi=sum((ATksi*covmat).*ATksi,2);
                        clear ATksi
                        
                        eta=(GM./(r.^2.*gamaP.*cos(fiG)).*sqrt(eta))*(180/pi)*3600;
                        eta(fi==pi/2 | fi==-pi/2)=0;
                        ksi=(GM./(r.^2.*gamaP).*sqrt(ksi))*(180/pi)*3600;

                        Pg=sqrt(eta.^2+ksi.^2);
                       
                        clear eta ksi
                    elseif volbapar(i)==5 %Disturbing potential                        
     
                        Pg=sum((AT*covmat).*AT,2);
                        clear AT
                        Pg=(GM./r.*sqrt(Pg));
                        
                    elseif volbapar(i)==6 %Disturbing tensor Trr_Tp_Tll
                        
                        Trr=sum((ATrr*covmat).*ATrr,2);
                        clear ATrr
                        Trr=(GM./r.^3.*sqrt(Trr))*10^9;
                        
                        Tfifi=sum((ATfifi*covmat).*ATfifi,2);
                        clear ATfifi
                        Tfifi=(GM./r.^3.*sqrt(Tfifi))*10^9;
                        
                        Tll=sum((ATll*covmat).*ATll,2);
                        clear ATll
                        Tll=(GM./(r.^3.*cos(fiG).^2).*sqrt(Tll))*10^9;
                        Tll(fi>pi/180*(89.9) | fi<pi/180*(-89.9))=0;
                        
                        Pg=Trr;
                        clear Trr
                        Pg=[Pg Tfifi];
                        clear Tfifi
                        Pg=[Pg Tll];
                        clear Tll
                        
                    elseif volbapar(i)==7 %Disturbing tensor Trp_Trl_Tpl
                        
                        Trfi=sum((ATrfi*covmat).*ATrfi,2);
                        clear ATrfi
                        Trfi=(GM./r.^3.*sqrt(Trfi))*10^9;
                        
                        Trl=sum((ATrl*covmat).*ATrl,2);
                        clear ATrl
                        Trl=(GM./(r.^3.*cos(fiG)).*sqrt(Trl))*10^9;
                        
                        Tfil=sum((ATfil*covmat).*ATfil,2);
                        clear ATfil
                        Tfil=(GM./(r.^3.*cos(fiG)).*sqrt(Tfil))*10^9;
                        Tfil(fi>pi/180*(89.5) | fi<pi/180*(-89.5),:)=0;
                        
                        Pg=Trfi;
                        clear Trfi
                        Pg=[Pg Trl];
                        clear Trl
                        Pg=[Pg Tfil];
                        clear Tfil
                        
                    elseif volbapar(i)==8 %Disturbing tensor Txx_Tyy_Tzz
                    elseif volbapar(i)==9 %Disturbing tensor Txy_Txz_Tyz
                    elseif volbapar(i)==10 %Geoid undulation
                        
                        Pg=sum((AN*covmat).*AN,2);
                        clear AN
                        Pg=(GM./(r.*gamaP).*sqrt(Pg));
                        
                    elseif volbapar(i)==11 %Gravitational potential
                        
                        Pg=sum((AV*covmat).*AV,2);
                        clear AV
                        Pg=(GM./r.*sqrt(Pg));
                        
                    elseif volbapar(i)==12 %Gravitational tensor Vrr_Vp_Vll
                        
                        Vrr=sum((AVrr*covmat).*AVrr,2);
                        clear AVrr
                        Vrr=(GM./r.^3.*sqrt(Vrr))*10^9;
                        
                        Vfifi=sum((AVfifi*covmat).*AVfifi,2);
                        clear AVfifi
                        Vfifi=(GM./r.^3.*sqrt(Vfifi))*10^9;
                        
                        Vll=sum((AVll*covmat).*AVll,2);
                        clear AVll
                        Vll=(GM./(r.^3.*cos(fiG).^2).*sqrt(Vll))*10^9;
                        Vll(fi>pi/180*(89.9) | fi<pi/180*(-89.9))=0;
                        
                        Pg=Vrr;
                        clear Vrr
                        Pg=[Pg Vfifi];
                        clear Vfifi
                        Pg=[Pg Vll];
                        clear Vll
                        
                    elseif volbapar(i)==13 %Gravitational tensor Vrp_Vrl_Vpl
                        
                        Vrfi=sum((AVrfi*covmat).*AVrfi,2);
                        clear AVrfi
                        Vrfi=(GM./r.^3.*sqrt(Vrfi))*10^9;
                        
                        Vrl=sum((AVrl*covmat).*AVrl,2);
                        clear AVrl
                        Vrl=(GM./(r.^3.*cos(fiG)).*sqrt(Vrl))*10^9;
                        
                        Vfil=sum((AVfil*covmat).*AVfil,2);
                        clear AVfil
                        Vfil=(GM./(r.^3.*cos(fiG)).*sqrt(Vfil))*10^9;
                        Vfil(fi>pi/180*(89.5) | fi<pi/180*(-89.5),:)=0;
                        
                        Pg=Vrfi;
                        clear Vrfi
                        Pg=[Pg Vrl];
                        clear Vrl
                        Pg=[Pg Vfil];
                        clear Vfil
                        
                    elseif volbapar(i)==14 %Gravitational tensor Vxx_Vyy_Vzz
                    elseif volbapar(i)==15 %Gravitational tensor Vxy_Vxz_Vyz
                    elseif volbapar(i)==16 %Gravity vector gX_gY_gZ
                        
                        Wr=sum((AWr*covmat).*AWr,2);
                        clear AWr
                        Wr=(GM./r.^2.*sqrt(Wr));
                        
                        Wfi=sum((AWfi*covmat).*AWfi,2);
                        clear AWfi
                        Wfi=(GM./r.^2.*sqrt(Wfi));
                        
                        Wlambda=sum((AWlambda*covmat).*AWlambda,2);
                        clear AWlambda
                        Wlambda=(GM./(r.^2.*cos(fiG)).*sqrt(Wlambda));
                                               
                        Pg=Wfi*10^5;
                        clear Wfi
                        Pg=[Pg Wlambda*10^5];
                        clear Wlambda
                        Pg=[Pg Wr*10^5];
                        clear Wr
                    elseif volbapar(i)==17 %Gravity sa
                        
                        Pg=sum((Ag_sa*covmat).*Ag_sa,2);
                        clear Ag_sa
                        Pg=(GM./r.^2.*sqrt(Pg))*10^5;
                        
                    elseif volbapar(i)==18 %Gravity potential
                        
                        Pg=sum((AW*covmat).*AW,2);
                        clear AW
                        Pg=(GM./r.*sqrt(Pg));
                        
                    elseif volbapar(i)==19 %Gravity anomaly sa
                        
                        Pg=sum((Aanomalia_sa*covmat).*Aanomalia_sa,2);
                        clear Aanomalia_sa
                        Pg=(GM./r.^2.*sqrt(Pg))*10^5;
                 
                    elseif volbapar(i)==20 %Gravity disturbance
                        
                        Wr_porucha=sum((AWr_porucha*covmat).*AWr_porucha,2);
                        clear AWr_porucha
                        Wr_porucha=(GM./r.^2.*sqrt(Wr_porucha));
                        
                        Wfi_porucha=sum((AWfi_porucha*covmat).*AWfi_porucha,2);
                        clear AWfi_porucha
                        Wfi_porucha=(GM./r.^2.*sqrt(Wfi_porucha));
                        
                        Wlambda_porucha=sum((AWlambda_porucha*covmat).*AWlambda_porucha,2);
                        clear AWlambda_porucha
                        Wlambda_porucha=(GM./(r.^2.*cos(fiG)).*sqrt(Wlambda_porucha));

                        Pg=sqrt(Wr_porucha.^2+Wlambda_porucha.^2+Wfi_porucha.^2)*10^5;

                        clear Wr_porucha Wfi_porucha Wlambda_porucha
                        
                    elseif volbapar(i)==21 %Gravity disturbance sa
                        
                        Pg=sum((Aporucha_sa*covmat).*Aporucha_sa,2);
                        clear Aporucha_sa
                        Pg=(GM./r.^2.*sqrt(Pg))*10^5;
                        
                    elseif volbapar(i)==22 %Height anomaly Ell
                        
                        Pg=sum((AzetaEl*covmat).*AzetaEl,2);
                        clear AzetaEl
                        Pg=(GM./(r.*gamaP).*sqrt(Pg));
     
                    elseif volbapar(i)==23 %Height anomaly
                        
                        Pg=sum((Azeta*covmat).*Azeta,2);
                        clear Azeta
                        Pg=(GM./(r.*gamaP).*sqrt(Pg));
                        
                    elseif volbapar(i)==24 %Second radial derivative of disturbing potential
                        
                        Pg=sum((AT_rr*covmat).*AT_rr,2);
                        clear AT_rr
                        Pg=(GM./r.^3.*sqrt(Pg))*10^9;
                        
                    elseif volbapar(i)==25 %Second radial derivative of gravity potential
                        
                        Pg=sum((AWrr*covmat).*AWrr,2);
                        clear AWrr
                        Pg=(GM./r.^3.*sqrt(Pg))*10^9;
                        
                    end

                    if i==1
                        P=Pg;
                        clear Pg
                    else
                        P=[P Pg];
                        if i==pocetpar
                           clear Pg
                        end
                    end                    
                end 
                
                clear q q2 gamaP covmat
            end
            
            set(findobj('tag','hlasky'),'string',...
                    '','fontsize',8,'foregroundcolor','k'); drawnow;

            cas=toc; %Stop clock

            %% Export data
            %==============================================================
                        
            if coord==1 %Spherical coordinates
                fi=fiG; %Ellipsoidal latitude is replaced by the spherical one
                clear fiG %Redundand vector fiG is deleted
            elseif coord==0 %Ellipsoidal coordinates
                clear fiG %Spherical latitude is deleted, since exported is
                %to be the ellipsoidal one
            end
            
            if Export_data_txt==1 || Export_data_mat==1 || nargout==1
                
                if GUI==1 %If working with the GUI
                    set(findobj('tag','hlasky'),'string',...
                        'Creating data file...','fontsize',8,...
                        'foregroundcolor','k'); drawnow;
                end
                
                if GUI==0 && Status_bar==1 %If working without the GUI
                    fprintf('Creating data file...\n');
                end
                
                if volbagridcheck==1
                    if STD==0
                        [lambda,fi]=meshgrid(180/pi*(lambda),180/pi*(fi));
                        [j_fi,i_fi]=size(fi);
                        [j_lambda,i_lambda]=size(lambda);
                        fi=fi(:);
                        lambda=lambda(:);
                    elseif STD==1
                        fi=180/pi*(fi);
                        lambda=180/pi*(lambda);
                    end
                elseif volbadiskcheck==1
                    fi=180/pi*(fi);
                    lambda=180/pi*(lambda);
                end                
                              
                try
                    if volbagridcheck==1
                        output=[fi';lambda';P'];
                    elseif volbadiskcheck==1
                        if coord==1 %Spherical coordinates
                            output=[fi';lambda';r';P'];
                        elseif coord==0 %Ellipsoidal coordinates
                            output=[fi';lambda';h';P'];
                        end
                    end

                    if Export_data_txt==1 %Export to txt file
                        [rows_output, cols_output]=size(output);

                        exp1=fopen([outadresar,[outname '.txt']],'w');
                        fprintf(exp1,...
                            [repmat('% 0.12e ',1,rows_output-1) '% 0.12e\n'],output);                   
                        fclose(exp1);
                        
                        set(findobj('tag','hlasky'),'string',...
                            '','fontsize',8,...
                            'foregroundcolor','k'); drawnow;  
                    end
    
                    if Export_data_mat==1 %Export to .mat file
                        output=output';
                        try
                            save([outadresar,[outname '.mat']],'output','-v7.3');
                        catch %#ok<CTCH>
                            save([outadresar,[outname '.mat']],'output');
                        end
                    end
                
                    if nargout==0
                        clear output
                    end
                    
                    if nargout==1 && Export_data_mat==0
                        output=output';
                    end
                    
                    export_error=0;
                catch err
                    export_error=1;
                    
                    try
                        save([outadresar,[outname '_phi.mat']],'fi','-v7.3');
                        save([outadresar,[outname '_lambda.mat']],'lambda','-v7.3');
                        save([outadresar,[outname '_functionals.mat']],'P','-v7.3');
                    catch %#ok<CTCH>
                        save([outadresar,[outname '_phi.mat']],'fi');
                        save([outadresar,[outname '_lambda.mat']],'lambda');
                        save([outadresar,[outname '_functionals.mat']],'P');
                    end
                end
            end
            
            %% Export report
            %============================================================== 
            
            if Export_report==1
                
                if GUI==1 %If working with the GUI
                    set(findobj('tag','hlasky'),'string',...
                        'Creating report file...','fontsize',8,...
                        'foregroundcolor','k'); drawnow;
                end
                
                if GUI==0 && Status_bar==1 %If working without the GUI
                    fprintf('Creating report file...\n')
                end
                
                date_time=clock;
                date_time=datestr(date_time);
                
                reportname=[outname '_Report' '.txt'];
                exp=fopen([outadresar,reportname],'w');
                                
                if Export_data_txt==0 && Export_data_mat==0
                    fi=180/pi*(fi);
                    lambda=180/pi*(lambda);
                    if STD==0
                        j_fi=length_fi;
                        i_fi=length(lambda);                    
                    end
                end  

                if GUI==1 %Working with the GUI
                    fprintf(exp,'Software                                           \t2.1.4 (GUI used)\n');
                elseif GUI==0 %Working without the GUI
                    fprintf(exp,'Software                                           \t2.1.4 (GUI not used)\n');
                end
                fprintf(exp,'Computer name                                      \t');
                if isunix
                    [ret,name]=system('hostname');
                    fprintf(exp,'%s\n',name);
                elseif ispc
                    fprintf(exp,'%s\n',getenv('computername'));
                else
                    fprintf(exp,'%s\n','');
                end
                fprintf(exp,'Generating date                                    \t');
                fprintf(exp,'%s',date_time(1:11));
                fprintf(exp,'\nGenerating time                                  \t');
                fprintf(exp,'%s',date_time(13:20));
                fprintf(exp,'\nMATLAB version                                   \t');
                fprintf(exp,'%s',version);
                fprintf(exp,'\nComputed                                         \t');
                
                if STD==0
                    fprintf(exp,'%s','Functionals of the geopotential');
                elseif STD==1
                    fprintf(exp,'%s','Commission errors');
                end
                
                if volbaloadcheck==1
                    fprintf(exp,'\nImported data file                                 \t');
                    fprintf(exp,'%s',[loadadresar,loadname]);
                end
                
                if STD==0
                    fprintf(exp,'\nGeopotential model file                          \t');
                    fprintf(exp,'%s',[GGMadresar,GGMname]);
                elseif STD==1
                    fprintf(exp,'\nVariance-covariance matrix file                  \t');
                    fprintf(exp,'%s',[GGMcovadresar,GGMcovname]);
                end
                
                if STD==0
                    if any(volbapar==10) || any(volbapar==23)
                        fprintf(exp,'\nDigital terrain model file                         \t');
                        fprintf(exp,'%s',[loadadresarDMR,loadnameDMR]);
                        fprintf(exp,'\nNewtonian gravitational constant (m^3*kg^-1*s^-2)  \t');
                        fprintf(exp,'%0.14e',G);
                        fprintf(exp,'\nDensity of the crust (kg*m^-3)                     \t');
                        fprintf(exp,'%0.14e',ro);
                    end
                end
                fprintf(exp,'\nGM of the geopotential model (m^3*s^-2)          \t');                
                fprintf(exp,'%0.14e',GM);
                fprintf(exp,'\nR of the geopotential model (m)                  \t');
                fprintf(exp,'%0.14e',R);
                fprintf(exp,'\nMinimum used degree                              \t');
                fprintf(exp,'%0.0f',nmin);
                fprintf(exp,'\nMaximum used degree                              \t');
                fprintf(exp,'%0.0f',nmax);  
                fprintf(exp,'\nReference ellipsoid                              \t');
                
                if length(ellipsoid)==1
                    if ellipsoid==1
                        fprintf(exp,'GRS80');
                    elseif ellipsoid==2
                        fprintf(exp,'WGS84');
                    end
                elseif length(ellipsoid)==5
                    fprintf(exp,'User defined (GM = %0.14e m^3*s^-2, a = %0.14e m , e = %0.14e, C20 = %0.14e, omega = %0.14e rad*sec^-1)',GMEl,aEl,eEl,CEl_20,omegaEl);
                end
                                
                fprintf(exp,'\nType of the input coordinates                    \t');
                if coord==1
                    fprintf(exp,'Spherical');
                elseif coord==0
                    fprintf(exp,'Ellipsoidal');
                end
                fprintf(exp,'\nLatitude limit North (deg)                       \t');
                fprintf(exp,'%0.14e',max(fi));
                fprintf(exp,'\nLatitude limit South (deg)                       \t');
                fprintf(exp,'%0.14e',min(fi));
                fprintf(exp,'\nLongitude limit West (deg)                       \t');
                fprintf(exp,'%0.14e',min(lambda));
                fprintf(exp,'\nLongitude limit East (deg)                       \t');
                fprintf(exp,'%0.14e',max(lambda));

                if volbagridcheck==1
                    fprintf(exp,'\nLatitude parallels                               \t');
                    fprintf(exp,'%0.0f',j_fi);
                    fprintf(exp,'\nLongitude parallels                              \t');
                    fprintf(exp,'%0.0f',i_fi);
                    fprintf(exp,'\nNumber of grid points                            \t');
                    fprintf(exp,'%0.0f',j_fi*i_fi);                                    
                    fprintf(exp,'\nGrid height above the reference surface (m)      \t');
                    
                    if coord==1 %Spherical coordinates
                        fprintf(exp,'%0.14e',hsph);
                    elseif coord==0 %Ellipsoidal coordinates
                        fprintf(exp,'%0.14e',h(1));
                    end
                    clear h
                    
                    hr=floor(cas/3600);
                    minutes=floor((cas/3600-hr)*60);
                    sec=round(((cas/3600-hr)*60-minutes)*60);
                    day_cpu=floor(hr/24);
                    hr=hr-day_cpu*24;
                    fprintf(exp,'\nComputation time (dd:hh:mm:ss)                     \t'); 
                    fprintf(exp,'%02d:%02d:%02d:%02d',day_cpu,hr,minutes,sec);
                    fprintf(exp,'\nComputation of fully normalized ALFs             \t');
                    
                    if volbaALFs==1
                        fprintf(exp,'Standard forward column method');
                    elseif volbaALFs==2
                        fprintf(exp,'Modified forward column method');
                    elseif volbaALFs==3
                        fprintf(exp,'Extended-range arithmetic');
                    end

                    fprintf(exp,'\n\nExported data file contains the following columns:');
                    
                    if coord==1 %Spherical coordinates
                        fprintf(exp,'\nSpherical Latitude (deg)   ');
                    elseif coord==0 %Ellipsoidal coordinates
                        fprintf(exp,'\nEllipsoidal Latitude (deg)   ');
                    end
                    fprintf(exp,'Longitude (deg)   ');                               

                elseif volbadiskcheck==1
                    fprintf(exp,'\nNumber of points                                 \t');
                    fprintf(exp,'%0.0f',length_fi);
                    hr=floor(cas/3600);
                    minutes=floor((cas/3600-hr)*60);
                    sec=round(((cas/3600-hr)*60-minutes)*60);
                    day_cpu=floor(hr/24);
                    hr=hr-day_cpu*24;
                    fprintf(exp,'\nComputation time (dd:hh:mm:ss)                     \t'); 
                    fprintf(exp,'%02d:%02d:%02d:%02d',day_cpu,hr,minutes,sec);
                    fprintf(exp,'\nComputation of fully normalized ALFs             \t');
                    
                    if volbaALFs==1
                        fprintf(exp,'Standard forward column method');
                    elseif volbaALFs==2
                        fprintf(exp,'Modified forward column method');
                    elseif volbaALFs==3
                        fprintf(exp,'Extended-range arithmetic');
                    end
                    
                    fprintf(exp,'\n\nExported data file contains the following columns:');
                    
                    if coord==1 %Spherical coordinates
                        fprintf(exp,'\nSpherical Latitude (deg)   ');
                    elseif coord==0 %Ellipsoidal coordinates
                        fprintf(exp,'\nEllipsoidal Latitude (deg)   ');
                    end
                    fprintf(exp,'Longitude (deg)   ');
                    
                    if coord==1
                        fprintf(exp,'Spherical radius (m)   ');
                    elseif coord==0
                        fprintf(exp,'Height above the reference ellipsoid (m)   ');
                    end
                    
                end
                               
                for i=1:pocetpar
                    if volbapar(i)~=1
                        if GUI==1 %If working with the GUI
                            Par=cellstr(get(findobj('tag',sprintf('P%0.0d',i)),'string'));
                        elseif GUI==0 %If working without the GUI
                            Par=eval(sprintf('P%0.0d',i));
                        end
                        Par=Par(volbapar(i));

                        fprintf(exp,'%s ',char(Par));
                        
                        if volbapar(i)==2 || volbapar(i)==3 
                            fprintf(exp,'(arcsec)   ');
                        elseif volbapar(i)==4
                            fprintf(exp,'(arcsec)   Azimuth (deg)   ');
                        elseif volbapar(i)==5 || volbapar(i)==11 || volbapar(i)==18
                            fprintf(exp,'(m^2*s^-2)   ');
                        elseif volbapar(i)==10 || volbapar(i)==22 || volbapar(i)==23
                            fprintf(exp,'(m)   ');
                        elseif volbapar(i)==16 || volbapar(i)==17 || volbapar(i)==19 || volbapar(i)==20 || volbapar(i)==21                        
                            fprintf(exp,'(mGal)   ');
                        elseif volbapar(i)==6 || volbapar(i)==7 || volbapar(i)==8 || volbapar(i)==9 || volbapar(i)==12 || volbapar(i)==13 || volbapar(i)==14 || volbapar(i)==15 || volbapar(i)==24 || volbapar(i)==25
                            fprintf(exp,'(Eotvos)   ');
                        end
                    end
                end
                
                if Export_data_txt==0 && Export_data_mat==0
                    export_error=0;
                end
                
                if export_error==1
                    fprintf(exp,'\n\nNote: GrafLab failed to create the data file due to the lack of memory.');
                    fprintf(exp,'\nHowever, GrafLab created three *.mat files, which contain all the computed data.');
                end
                
                fclose(exp); 
                
                set(findobj('tag','hlasky'),'string',...
                    '','fontsize',8,'foregroundcolor','k'); drawnow; 
            end
            
            %% Display data
            %==============================================================
            
            if (display_data==1 || display_data==2) && volbagridcheck==1                                
                
                if STD==0
                    if Export_data_txt==1 || Export_data_mat==1
                        fiGrid=reshape(fi,length_fi,length_lambda);
                        lambdaGrid=reshape(lambda,length_fi,length_lambda);
                    elseif Export_report==1
                        if nargout==0
                            [lambdaGrid,fiGrid]=meshgrid(lambda,fi);
                        else
                            fiGrid=reshape(fi,length_fi,length_lambda)*pi/180;
                            lambdaGrid=reshape(lambda,length_fi,length_lambda)*pi/180;
                        end
                    else
                        if nargout==0
                            fi=180/pi*(fi);
                            lambda=180/pi*(lambda);
                            [lambdaGrid,fiGrid]=meshgrid(lambda,fi);
                        else
                            fiGrid=reshape(fi,length_fi,length_lambda);
                            lambdaGrid=reshape(lambda,length_fi,length_lambda);
                        end
                    end
                elseif STD==1
                    length_fi=length_fi_plot;
                    if Export_data_txt==1 || Export_data_mat==1
                        fiGrid=reshape(fi,length_fi,length_lambda);
                        lambdaGrid=reshape(lambda,length_fi,length_lambda);
                    elseif Export_report==1
                        fiGrid=reshape(fi,length_fi,length_lambda)*pi/180;
                        lambdaGrid=reshape(lambda,length_fi,length_lambda)*pi/180;
                    else
                        fiGrid=reshape(fi,length_fi,length_lambda);
                        lambdaGrid=reshape(lambda,length_fi,length_lambda);
                    end
                end
                
                clear fi lambda
                
                %Colormap
                if colmap==1
                    colmap='jet';
                elseif colmap==2
                    colmap='hsv';
                elseif colmap==3
                    colmap='hot';
                elseif colmap==4
                    colmap='cool';
                elseif colmap==5
                    colmap='spring';
                elseif colmap==6
                    colmap='summer';
                elseif colmap==7
                    colmap='autumn';
                elseif colmap==8
                    colmap='winter';
                elseif colmap==9
                    colmap='gray';
                elseif colmap==10
                    colmap='bone';
                elseif colmap==11
                    colmap='copper';
                elseif colmap==12
                    colmap='pink';
                elseif colmap==13
                    colmap='lines';
                end
                
                %Graphic format file
                if volbaformat==1
                    format='bmp';
                    dformat='-dbmp16m';
                elseif volbaformat==2
                    format='emf';
                    dformat='-dmeta';
                elseif volbaformat==3
                    format='eps';
                    dformat='-depsc';
                elseif volbaformat==4
                    format='jpeg';
                    dformat='-djpeg';
                elseif volbaformat==5
                    format='pdf';
                    dformat='-dpdf';
                elseif volbaformat==6
                    format='png';
                    dformat='-dpng';
                elseif volbaformat==7
                    format='tiff';
                    dformat='-dtiff';
                end
                
                if display_data==1
                    coast=load('coast.mat'); %Loading continents  
                end
                
                zobr=1;
      
                for i=1:pocetpar
                    
                    if GUI==1 %If working with the GUI
                        if i==1
                            set(findobj('tag','hlasky'),'string',...
                                'Displaying 1st functional...','fontsize',...
                                8,'foregroundcolor','k'); drawnow;
                        elseif i==2
                            set(findobj('tag','hlasky'),'string',...
                                'Displaying 2nd functional...','fontsize',...
                                8,'foregroundcolor','k'); drawnow;
                        elseif i==3
                            set(findobj('tag','hlasky'),'string',...
                                'Displaying 3rd functional...','fontsize',...
                                8,'foregroundcolor','k'); drawnow;
                        elseif i==4
                            set(findobj('tag','hlasky'),'string',...
                                'Displaying 4th functional...','fontsize',...
                                8,'foregroundcolor','k'); drawnow;
                        end
                    end
                    
                    if GUI==0 && Status_bar==1 %If working without the GUI
                        if i==1
                            fprintf('Displaying 1st functional...\n')
                        elseif i==2
                            fprintf('Displaying 2nd functional...\n')
                        elseif i==3
                            fprintf('Displaying 3rd functional...\n')
                        elseif i==4
                            fprintf('Displaying 4th functional...\n')
                        end
                    end
                                                                              
                    if volbapar(i)~=1                                     
                        if volbapar(i)==6 || volbapar(i)==7 || volbapar(i)==8 || volbapar(i)==9 || volbapar(i)==12 || volbapar(i)==13 || volbapar(i)==14 || volbapar(i)==15 || volbapar(i)==16
                            for j=0:2
                                Pdisp=reshape(P(:,zobr+j),length_fi,length_lambda);

                                if volbapar(i)==6
                                    nazov=['Disturbing_tensor_Trr';...
                                        'Disturbing_tensor_Tpp';...
                                        'Disturbing_tensor_Tll'];
                                elseif volbapar(i)==7
                                    nazov=['Disturbing_tensor_Trp';...
                                        'Disturbing_tensor_Trl';...
                                        'Disturbing_tensor_Tpl'];
                                elseif volbapar(i)==8
                                    nazov=['Disturbing_tensor_Txx';...
                                        'Disturbing_tensor_Tyy';...
                                        'Disturbing_tensor_Tzz'];
                                elseif volbapar(i)==9
                                    nazov=['Disturbing_tensor_Txy';...
                                        'Disturbing_tensor_Txz';...
                                        'Disturbing_tensor_Tyz'];
                                elseif volbapar(i)==12
                                    nazov=['Gravitational_tensor_Vrr';...
                                        'Gravitational_tensor_Vpp';...
                                        'Gravitational_tensor_Vll'];
                                elseif volbapar(i)==13
                                    nazov=['Gravitational_tensor_Vrp';...
                                        'Gravitational_tensor_Vrl';...
                                        'Gravitational_tensor_Vpl'];
                                elseif volbapar(i)==14
                                    nazov=['Gravitational_tensor_Vxx';...
                                        'Gravitational_tensor_Vyy';...
                                        'Gravitational_tensor_Vzz'];
                                elseif volbapar(i)==15
                                    nazov=['Gravitational_tensor_Vxy';...
                                        'Gravitational_tensor_Vxz';...
                                        'Gravitational_tensor_Vyz'];
                                elseif volbapar(i)==16
                                    nazov=['Gravity_vector_gX';...
                                        'Gravity_vector_gY';...
                                        'Gravity_vector_gZ'];
                                end

                                Okno=figure('numbertitle','off','name',...
                                    char(nazov(j+1,:)),'visible','off');
                                
                                if display_data==1 %Uses Mapping Toolbox
                                    %Note that the
                                    %processing time may take a long time in
                                    %case of large data sets.
                                    worldmap([min(min(fiGrid)) max(max(fiGrid))],...
                                        [min(min(lambdaGrid)) max(max(lambdaGrid))]);
                                    contourfm(fiGrid,lambdaGrid,Pdisp,ncolor,...
                                        'linestyle','none');
                                elseif display_data==2 %Does not use Mapping Toolbox
                                    %A significantly faster way how to display
                                    %large data sets is to use the command
                                    %"imagesc". However, note that this
                                    %function displays the grid data simply as
                                    %a matrix, that is, common cartographic
                                    %projections are not supported here.
                                    imagesc(lambdaGrid(1,:),fiGrid(:,1),Pdisp)
                                    set(gca,'Ydir','normal')
                                    xlabel('Longitude (deg)');
                                    ylabel('Latitude (deg)');
                                    grid on;
                                end
                                colormap(sprintf('%s(%d)',colmap,ncolor));
                                labelcolbar=colorbar('location','southoutside');
                                mean_Pdisp=mean(mean(Pdisp));
                                std_Pdisp=std(Pdisp(:));
                                caxis([max([mean_Pdisp-4*std_Pdisp min(min(Pdisp))]) min([mean_Pdisp+4*std_Pdisp max(max(Pdisp))])]);

                                %Units in colorbar
                                if volbapar(i)==16
                                    set(get(labelcolbar,'xlabel'),'string','mGal');
                                else
                                    set(get(labelcolbar,'xlabel'),'string','Eotvos');
                                end
                                
                                if display_data==1
                                    geoshow(coast.lat,coast.long,'color','black');
                                end
                                
                                if isempty(outadresar)
                                    print(Okno,sprintf('%s',dformat),...
                                        sprintf('-r%i',DPI),...
                                        sprintf('%s_%s.%s',outname,char(nazov(j+1,:)),format));
                                else
                                    print(Okno,sprintf('%s',dformat),...
                                        sprintf('-r%i',DPI),...
                                        sprintf('%s%s%s_%s.%s',outadresar,filesep,...
                                        outname,char(nazov(j+1,:)),format));
                                end
                            end                                                       
                            
                            zobr=zobr+3;                          
                        else                                
                            Pdisp=reshape(P(:,zobr),length_fi,length_lambda);
                            
                            %Deleting azimuth if Theta has been computed
                            if volbapar(i)==4 && STD==0
                                P(:,zobr+1)=[];
                            end
                            
                            if GUI==1 %If working with the GUI
                                Nazov_okna=cellstr(get(findobj('tag','P1'),'string'));
                            elseif GUI==0 %If working without the GUI
                                Nazov_okna=P1;
                            end
                            Nazov_okna=Nazov_okna(volbapar(i));                            
                            Okno=figure('numbertitle','off','name',...
                                char(Nazov_okna),'visible','off');
                            
                            if display_data==1 %Uses Mapping Toolbox
                                %Note that the
                                %processing time may take a long time in
                                %case of large data sets.
                                worldmap([min(min(fiGrid)) max(max(fiGrid))],...
                                    [min(min(lambdaGrid)) max(max(lambdaGrid))]);
                                contourfm(fiGrid,lambdaGrid,Pdisp,ncolor,...
                                    'linestyle','none');
                            elseif display_data==2 %Does not use Mapping Toolbox
                                %A significantly faster way how to display
                                %large data sets is to use the command
                                %"imagesc". However, note that this
                                %function displays the grid data simply as
                                %a matrix, that is, common cartographic
                                %projections are not supported here.
                                imagesc(lambdaGrid(1,:),fiGrid(:,1),Pdisp)
                                set(gca,'Ydir','normal')
                                xlabel('Longitude (deg)');
                                ylabel('Latitude (deg)');
                                grid on;
                            end
                            
                            colormap(sprintf('%s(%d)',colmap,ncolor));
                            labelcolbar=colorbar('location','southoutside');
                            mean_Pdisp=mean(mean(Pdisp));
                            std_Pdisp=std(Pdisp(:));
                            caxis([max([mean_Pdisp-4*std_Pdisp min(min(Pdisp))]) min([mean_Pdisp+4*std_Pdisp max(max(Pdisp))])]);
                          
                            %Units in colorbar
                            if volbapar(i)==2 || volbapar(i)==3 || volbapar(i)==4
                                set(get(labelcolbar,'xlabel'),'string','arcsec');
                            elseif volbapar(i)==5 || volbapar(i)==11 || volbapar(i)==18
                                set(get(labelcolbar,'xlabel'),'string','m^2 \cdot s^{-2}');
                            elseif volbapar(i)==10 || volbapar(i)==22 || volbapar(i)==23
                                set(get(labelcolbar,'xlabel'),'string','m');
                            elseif volbapar(i)==16 || volbapar(i)==17 || volbapar(i)==19 || volbapar(i)==20 || volbapar(i)==21                        
                                set(get(labelcolbar,'xlabel'),'string','mGal');
                            elseif volbapar(i)==6 || volbapar(i)==7 || volbapar(i)==8 || volbapar(i)==9 || volbapar(i)==12 || volbapar(i)==13 || volbapar(i)==14 || volbapar(i)==15 || volbapar(i)==24 || volbapar(i)==25
                                set(get(labelcolbar,'xlabel'),'string','Eotvos');
                            end

                            if display_data==1
                                geoshow(coast.lat,coast.long,'color','black');   
                            end
                            
                            if isempty(outadresar)
                                print(Okno,sprintf('%s',dformat),...
                                    sprintf('-r%i',DPI),sprintf('%s_%s.%s',...
                                    outname,char(Nazov_okna),format));
                            else
                                print(Okno,sprintf('%s',dformat),...
                                    sprintf('-r%i',DPI),sprintf('%s%s%s_%s.%s',...
                                    outadresar,filesep,outname,char(Nazov_okna),format));
                            end
                            
                            zobr=zobr+1;
                        end  
                    end
                    close(Okno)
                end
    
                set(findobj('tag','hlasky'),'string',...
                    '','fontsize',8,'foregroundcolor','k'); drawnow;
            end 
            
            set(findobj('tag','R_text'),'userdata','');
            set(findobj('tag','ell_text'),'userdata','');
          
            set(findobj('tag','OutFile'),'string','');
            
            if GUI==1 %Working with the GUI
                set(findobj('tag','hlasky'),'string',...
                    'Computation has been finished','fontsize',8,...
                    'foregroundcolor','k'); drawnow; 
            end
            
            if GUI==0 && Status_bar==1 %Working without the GUI
                fprintf('Computation has been finished\n')
            end
                
        case('Close')
            close all
    end
end

if nargout
    varargout{1}=output;
end

%SUBFUNCTIONS
%==========================================================================
%--------------------------------------------------------------------------
function [X,Y,Z]=sph2cart(lambda,fi,r)
%This function transforms spherical coordinates into cartesian coordinates.
%
%INPUT: "lambda" -- vector of spherical longitudes in radians
%       "fi"     -- vector of spherical latitudes in radians
%       "r"      -- vector spherical radii in metres
%
%OUTPUT: "X","Y","Z" -- three vectors of cartesian coordinates in metres

X=r.*cos(fi).*cos(lambda);
Y=r.*cos(fi).*sin(lambda);
Z=r.*sin(fi);
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
function [lambda,fi,r]=cart2sph(X,Y,Z) %#ok<DEFNU>
%This function transforms cartesian coordinates into spherical coordinates.
%
%INPUT: "X","Y","Z" -- three vectors of cartesian coordinates in metres
%
%OUTPUT: "lambda" -- vector of spherical longitudes in radians
%        "fi"     -- vector of spherical latitudes in radians
%        "r"      -- vector spherical radii in metres

r=sqrt(X.^2+Y.^2+Z.^2);
fi=atan(Z./sqrt(X.^2+Y.^2));
lambda=atan2(Y,X);
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
function [X,Y,Z]=ell2cart(fi,lambda,h,ellipsoid)
%This function transforms ellipsoidal coordinates into cartesian coordinates.
%
%INPUT: "fi"        -- vector of ellipsoidal latitudes in radians
%       "lambda"    -- vector of ellipsoidal longitudes in radians
%       "h"         -- vector ellipsoidal heights in metres
%       "ellipsoid" -- vector defining the semimajor axis of the
%                      reference ellipsoid and its first numerical 
%                      eccentricity, respectivelly.
%
%OUTPUT: "X","Y","Z" -- three vectors of cartesian coordinates

aEl=ellipsoid(1);
eEl=ellipsoid(2);

N=aEl./sqrt(1-eEl^2.*sin(fi).^2);

X=(N+h).*cos(fi).*cos(lambda);
Y=(N+h).*cos(fi).*sin(lambda);
Z=((1-eEl^2).*N+h).*sin(fi);
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
function [fi,lambda,h]=cart2ell(X,Y,Z,ellipsoid)
%This function transforms cartesian coordinates into ellipsoidal coordinates.
%
%INPUT: "X","Y","Z" -- three vectors of cartesian coordinates in metres
%       "ellipsoid" -- vector defining the semimajor axis of the
%                      reference ellipsoid and its first numerical 
%                      eccentricity, respectively.
% 
%OUTPUT: "fi"     -- vector of ellipsoidal latitudes in radians
%        "lambda" -- vector of ellipsoidal longitudes in radians
%        "h"      -- vector ellipsoidal heights in metres

aEl=ellipsoid(1);
eEl=ellipsoid(2);

lambda=atan2(Y,X);

p=sqrt(X.^2+Y.^2);
fi_temp=atan(Z./((1-eEl^2).*p));
lat_diff=9999;

while lat_diff>eps
    
    N=aEl./sqrt(1-eEl^2.*sin(fi_temp).^2);
    h=p./cos(fi_temp)-N;
    fi=atan(Z./((1-eEl^2.*N./(N+h)).*p));
    
    lat_diff=max(abs(fi-fi_temp));
    
    fi_temp=fi;

end

%--------------------------------------------------------------------------
%==========================================================================
