clear all ;
clc ;
% Code uses the TMM method t o c a l c u l a t e r e f l e c t i o n and transm iss ion through
% a s t r u c t u r e de f ined by TMMstructure using the propagat ion and s c a t t e r i n g
% matrix o f each per iod/ i n t e r f a c e .
StartLambda = 1000 ; EndLambda = 1500 ; L_diff = 0.1 ;
wavelength = StartLambda : L_diff : EndLambda ;
dispersion = 0 ; % ma ter ial  p r o p e r t i e s    f i x e d or wavelength dependent
MyMaterialBuilder(dispersion, wavelength ) ;
    % Design parameters
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    N = 18 ; % NumberOfPeriods in DBR (1/2 more than t o t a l   i . e . 18 −> 17.5 pa irs )
    MetalDepth = 25 ;
    Spacer1Depth = 0 ;
    Spacer2Depth = 75 ;
    NoOfLayers = 2*N + 6 ;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    TE_polarized = 1 ; % TE p o la r iz ed == 1 , TM p o la r iz ed == 0
    phi_0 = asin (0.0) ; % radians
    reverse = 0 ; % excitationside ( metal first == 0 , substratefirst == 1 )

 % The procedure i s adapted from J . Appl . Phys Vol 86 , No. 1 (1999 ) p.487
% and JAP 93 No. 7 p . 3693.
plotLambda = 1300 ;
if isempty(find( plotLambda == wavelength , 1 ) )
    error ( 'p l o t t i n g wavelength out o f range ' )
end
l = find ( plotLambda == wavelength ) ;
% Build s t r u c t u r e
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if dispersion == 0 % no d i sp e r s i on ;
    TMMstructure = BuildStructure (N, MetalDepth , Spacer1Depth , Spacer2Depth , 1 , reverse ) ;
else
    TMMstructure = BuildStructure (N, MetalDepth , Spacer1Depth , Spacer2Depth , l , reverse ) ;
end
t = zeros ( 1 , NoOfLayers ) ;
for layer = 1: NoOfLayers
    t (layer) = TMMstructure(layer) . depth ;
end
t_sum = cumsum( t ) ;

z = 0.5 : 1 : t_sum ( NoOfLayers ) ;
z_mat = sum( repmat ( z , NoOfLayers , 1 ) > repmat ( t_sum'  , 1 , length ( z ) ) , 1 ) + 1 ;
n = zeros ( 1 ,length ( z ) ) ;
for z_i = 1:length ( z )
    n ( z_i ) = TMMstructure ( z_mat ( z_i ) ) . ncomplex ;
end


    for l = 1:length ( wavelength )
        if dispersion == 1
            TMMstructure = BuildStructure (N, MetalDepth , Spacer1Depth , Spacer2Depth , l , reverse ) ;
        end
        n01 = TMMstructure( 1 ).ncomplex ; TMMstructure ( 1 ).q = cos ( phi_0 ) ;
        n02 = TMMstructure( NoOfLayers ).ncomplex ; TMMstructure ( NoOfLayers ).q = sqrt ( 1-( n01*sin ( phi_0 ) / n02 )^ 2 ) ;
        % For each layer , ob ta in   the   l a y e r   propagat ion   matrix ( L_mat )
        k0 = (2*pi / wavelength ( l ) ) ;
        for layer = 2: ( NoOfLayers-1)
            TMMstructure ( layer ) . q = sqrt ( 1- ( n01*sin ( phi_0 ) / TMMstructure ( layer ) . ncomplex )^ 2 ) ; % c os ( ph i _ la yer )
            xi = k0 * TMMstructure ( layer ) . ncomplex * TMMstructure ( layer ) . q ; % e f f e c t i v e wavevec tor
            TMMstructure ( layer ) . L_mat = [ exp(-1i * xi *TMMstructure ( layer ) . depth ) 0 ; 0 exp(1i * xi *TMMstructure ( layer ) . depth ) ] ;
        end

        % For each i n t e r f a c e , c a l c u l a t e the i n t e r f a c e matrix ( I_mat )
        NoOfInterfaces = NoOfLayers-1 ;
        for layer = 1: NoOfInterfaces % i n t e r f a c e ( j ) i s the i n t e r f a c e between l a y e r j and l a y e r j +1
            n1 = TMMstructure ( layer ) . ncomplex ; n2 = TMMstructure ( layer + 1 ). ncomplex ;
            q1 = TMMstructure ( layer ) . q ; q2 = TMMstructure ( layer + 1 ). q ;
            if TE_polarized == 1 % s e e SPIE 7521 , 75210G ( 2 0 0 9 ).
                interfaces(layer) . r_mat = ( n1*q1 - n2*q2 ) / ( n1*q1 + n2*q2 ) ;
                interfaces(layer) . t_mat = (2*n1*q1 ) / ( n1*q1 + n2*q2 ) ;
            elseif TE_polarized == 0
                interfaces(layer) . r_mat = ( n2*q1 - n1*q2 ) / ( n2*q1 + n1*q2 ) ;
                interfaces(layer) . t_mat = (2*n1*q1 ) / ( n2*q1 + n1*q2 ) ;
            end
           interfaces(layer) . I_mat = ( 1 /interfaces(layer) . t_mat )* [1  interfaces(layer) . r_mat ;  interfaces(layer) . r_mat  1 ] ;
        end
        
        S_0 = interfaces( NoOfInterfaces ).I_mat ;
        for m = NoOfInterfaces : -1:2
            S_0 = interfaces (m-1 ).I_mat * TMMstructure (m) . L_mat * S_0 ;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % From the f i n a l S matrix , c a l c u l a t e the r e f l e c t i o n and the
        % transm iss ion c o e f f i c i e n t s ( f o r t h i s wavelength , l )
        % r e f l e c t i o n R=|r|^2 = r * c on j ( r )
        r = S_0 ( 2 , 1 ) / S_0 ( 1 , 1 ) ;
        reflection( l ) = r*conj( r ) ;
        % transm iss ion T=| t |^2 * ( n _subs tra te / n_ambient_medium )
        % d e f in e transm iss i on _sca le i s r a t i o o f r e f r a c t i v e i n d i c i e s n _subs tra te / n_ambient_medium
        % * [ c os ( ph i _ sub s t ra t e ) / c os ( phi_0 ) ]
        transmission_ratio = ( n02 * TMMstructure ( NoOfLayers ) . q ) / ( n01 * TMMstructure ( 1 ) . q ) ;
        transmission( l ) = abs ( 1 / S_0 ( 1 , 1 ) ) ^ 2 * real ( transmission_ratio) ;
        %abs orp t i on ( l ) = 1 − ( r e f l e c t i o n ( l ) + transm iss ion ( l ) ) ;
        % Ca lcu la te the r e f l e c t i o n phase
        S_1 =interfaces( NoOfInterfaces-1).I_mat ;
        for m = ( NoOfInterfaces -1 ): -1:3
            S_1 = interfaces (m-1 ). I_mat * TMMstructure (m) . L_mat * S_1 ;
        end
        phaseangle( l )= angle( S_1 ( 2 , 1 ) / S_1 ( 1 , 1 ) ) ;
        if mod ( l , ( 100 / L_diff ) ) == 1
            disp ( [ 'Lambda = '  num2str( wavelength ( l ) )  'nm' ] )
        end
    end 

for layer = 2: NoOfLayers
    S_prime = interfaces (layer -1). I_mat ;
    for m = ( layer -1 ): -1:2
        S_prime =interfaces (m-1 ). I_mat * TMMstructure (m) . L_mat * S_prime ;
    end
    S_primeprime = interfaces ( NoOfInterfaces ) . I_mat ;
    for m = NoOfInterfaces : -1: ( layer +1)
        S_primeprime = interfaces(m-1 ). I_mat * TMMstructure (m) . L_mat * S_primeprime ;
    end
    z_f = find ( z_mat == layer ) ;
    % in d i c e s o f z array tha t are within tha t p a r t i cu l a r l a y e r
    z_pos = z ( z_f ) - t_sum ( layer -1) ;
    % p o s i t i o n withn l a y e r ( r e l a t i v e t o ( la yer −1) i n t e r f a c e
    sj = t ( layer ) ; sj_prime = sj - z_pos ;
    zeta_j = 2*pi / plotLambda * TMMstructure ( layer ) . ncomplex * TMMstructure ( layer ) . q ;
    E( z_f ) = ( S_primeprime ( 1 , 1 ) *exp(-1i * zeta_j * sj_prime ) + S_primeprime ( 2 , 1 ) *exp(1i * zeta_j * sj_prime ) )./ ...
     ( S_prime ( 1 , 1 ) * S_primeprime ( 1 , 1 ) *exp(-1i * zeta_j * sj ) + S_prime ( 1 , 2 ) * S_primeprime ( 2 , 1 ) *exp(1i * zeta_j * sj ) ) ;
end

E_p = abs (E )'  ;
plot ( z , n , z , abs (E ) )




    
function MyMaterialBuilder ( dispersion , wavelength )
% makes ma ter ia l f i l e s using the data g iven in the t a b l e s in F i l eD i r e c t o r y
% taken from r e f e r e n c e s ( e . g . Johnson and Chr is t y )
% check f i l e d i r e c t o r y in forma t ion i s c o r r e c t
FileDirectory = 'Z:\ ' ;
FilePath = 'Z:\ ' ;
% check names o f ma ter ia l f i l e s match
n_1 = 'GaAs ' ; n_2 = ' AlAs ' ;
n_metal = ' gold ' ;
f_1 = fopen([FilePath n_1 ' . t x t '],'wt') ;
f_2 = fopen ([FilePath n_2 ' . t x t '],'wt') ;
f_metal = fopen ( [ FilePath n_metal ' . t x t ' ] , 'wt' ) ;
if dispersion == 0
    fprintf ( f_1 , '%f %f %f ' , 1 , 3.5 , 0 ) ;
    fprintf ( f_2 , '%f %f %f ' , 1 , 2.9 , 0 ) ;
    fprintf ( f_metal , '%f %f %f ' , 1 , 0.04 , 6.312 ); % r e a l and imaginery component
else
    f_1 = fopen ( [ FilePath n_1 ' . t x t ' ] ,'wt') ;
    f_2 = fopen ( [ FilePath n_2 ' . t x t ' ] , 'wt'  ) ;
    f_metal = fopen ( [FilePath n_metal ' . t x t ' ] , 'wt' ) ;
    %{
    % i n t e r p o l a t i o n −− need t o i n t e r p o l a t e r e f r a c t i v e in d i c e s so they are
    % un i form ly sampled as the same frequenc y o f the chosen wavelength using
    % in te rp1 ( x , v ( x ) , qp )
    % x , the sample p o in t s ( wavelength )
    % v ( x ) , the va lues wanting t o be i n t e r p o l a t e d ( r e f r a c t i v e index )
    % qp , query p o in t s − the uniform spac ing wanted f o r x ( wavelength )
    %}
    A = load ( [ FileDirectory n_1 ' . t x t' ] ) ;
    n_GaAs = interp1(1000*A ( : , 1 ) ,A ( : , 2 ) , wavelength , 'spline' ) + ...
        1i *interp1(1000*A ( : , 1 ) ,A ( : , 3 ) , wavelength , 'spline') ;
    %n_GaAs = in t erp1 (1000*A ( : , 1 ) , A ( : , 2 ) , wavelength , ’ sp l ine ’ ) ;
    fprintf ( f_1 , '%f %f %f \n ' , [ wavelength ; real ( n_GaAs ) ; imag( n_GaAs ) ] ) ;
    A = load ( [FileDirectory n_2 ' . t x t ' ] ) ; % c on ta ins r e a l par ts on ly
    n_AlAs = interp1(1000*A ( : , 1 ) ,A ( : , 2 ) , wavelength , ' s p l in e ' ) ;
    fprintf ( f_2 , '%f %f \n ' , [ wavelength ; n_AlAs ] ) ;
    A = load ( [ FileDirectory n_metal ' . t x t ' ] ) ;
    n_metal = interp1(1000*A ( : , 1 ) ,A ( : , 2 ) , wavelength , ' s p l in e ' ) + ...
    1i*interp1(1000*A ( : , 1 ) ,A ( : , 3 ) , wavelength , ' s p l in e ' ) ;
    %n_metal = 0 + 1 i * in te rp1 (1000*A ( : , 1 ) , A ( : , 3 ) , wavelength , ’ sp l ine ’ ) ;
    fprintf ( f_metal , ' %f %f %f \n ' , [ wavelength ; real ( n_metal ) ; imag( n_metal ) ] ) ;
end
fclose ( 'all' ) ;
end

function TMMstructure = BuildStructure ( NoOfPeriods , MetalDepth , Spacer1Depth , Spacer2Depth , lambda , reverse )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% c r e a t e s a s t r u c t data t ype ( ’ TMMstructure ’ ) c on ta in ing the s t r u c t u r e t o be simulated .
% I t c on ta ins Q e n t r i e s , where Q i s the number o f la yers , and (Q−1) i s the number o f i n t e r f a c e s
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% s e l e c t s t r u c t u r e ma ter ia ls here . Make sure tha t a ) Ma ter ia lBu i lder has
% been updated and tha t ma ter ia l f i l e s e x i s t s in \ma ter ia ls f o l d e r b ) f i l e
% d i r e c t o r y ( e . g . F: \ ) i s c o r r e c t
FileDirectory = 'Z:\ ' ;
N_1 = 'GaAs ' ; N_2 = ' AlAs ' ;
N_metal = ' gold ' ;

A = load ( [ FileDirectory N_1 ' . t x t ' ] ) ;
n_1 = A( lambda , 2 ) + 1i *A( lambda , 3 ) ;

A = load ( [FileDirectory N_2 ' . t x t ' ] ) ;
n_2 = A( lambda , 2 ) ;

A = load ( [ FileDirectory N_metal ' . t x t ' ] ) ;
n_metal = A( lambda , 2 ) + 1i *A( lambda , 3 ) ;

% The semi−i n f i n i t e input and su b s t r a t e l a y e r ( normally a i r and GaAs )
n_air = 1 ; n_substrate = n_1 ;
n_spacer = 1.4 ;

Q = 2*NoOfPeriods + 6 ; % To ta l number o f l a y e r s

% t h i c k n e s s e s (nm) o f DBR pa ir ( [ h i gh e r _ r e f l ow e r _ r e f ] e . g . [GaAs AlAs ] ) )
H = [95 110 ];
% l a y e r q = 1 : the semi−i n f i n i t e a i r l a y e r
TMMstructure ( 1 ) . depth = 0; % depth , a r e a l number r e p r e s en t in g the th i c kn e s s o f n^th l a y e r in nm
TMMstructure ( 1 ) . ncomplex = n_air ; % ncomplex , a complex number correspond ing t o the r e f r a c t i v e index

TMMstructure ( 2 ) . depth = 1000 ;
TMMstructure ( 2 ) . ncomplex = n_air ;

% l a y e r q = 3 : the metal l a y e r
TMMstructure ( 3 ) . depth = MetalDepth ;
TMMstructure ( 3 ) . ncomplex = n_metal ;

%l a y e r q = 4 : the spacer l a y e r
TMMstructure ( 4 ) . depth = Spacer1Depth ;
TMMstructure ( 4 ) . ncomplex = n_spacer ;

TMMstructure ( 5 ) . depth = Spacer2Depth ;
%TMMstructure ( 5 ) . depth = H( 1 ) ;
TMMstructure ( 5 ) . ncomplex = n_1 ;
% l a y e r 6 −> Q−2 : the Bragg mu l t i la yer
for q = 6 : (Q-2)
    if mod ( q , 2 ) == 1
        TMMstructure ( q ) . depth = H( 1 ) ;
        TMMstructure ( q ) . ncomplex = n_1 ;
    else
        TMMstructure ( q ) . depth = H( 2 ) ;
        TMMstructure ( q ) . ncomplex = n_2 ;
    end
end

TMMstructure (Q-1 ). ncomplex = n_substrate ;
TMMstructure (Q-1 ). depth = 1000;

TMMstructure (Q ) . ncomplex = n_substrate ;
TMMstructure (Q ) . depth = 0;

if reverse == 1 % r e v e r s e s order o f s t r u c t u r e ( i . e . e x c i t e s from su b s t r a t e s id e )
    TMMstructure2 = TMMstructure ;
    for j = 1:Q
        TMMstructure ( j ) . ncomplex = TMMstructure2 (Q-j + 1 ). ncomplex ;
        TMMstructure ( j ) . depth = TMMstructure2 (Q-j + 1 ). depth ;
    end
    clear TMMstructure2
end

end




