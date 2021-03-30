"""
    expected_number_of_cluster_2PD(n::N, θ::T, σ::T, ntrunc::N)

Compute the expected number of clusters for a 2-parameter Poisson Dirichlet prior, by direct computation up to ntrunc. If ntrunc = n, this is an exact computation, if it turns out to be too long and that you expect large values of k to have a very small mass, you can set ntrunc < n.

# Examples
```julia-repl
julia> GibbsTypePriors.expected_number_of_cluster_2PD(100, 1.0, 0.1, 50)
[6.66854572439020053533515735851452973224981540358173684266962825487333967458307765332656984810574842651474015538338720078320799873773484654474045662876311930900255401608438351708575902470003695195472950846209930364916574756258195906978792735380469732399701637049873218563662828748041941268348464936340484967657964701446269020409816292488646213414168246181455914819473889600572869116478059869754289011725089446116537699336280040328071818457251209457817657352446105060779361083211752873341183127245705781828124887118172387163760512916480460420221393531464292989217584922021599146160085328610112860589003070634636891323889965306715193866013221311176472638005617254383323549220868258800194958774473369822418060972802272418739639187815136375304624761993600634276521918590366043000203183869628505847657460467697757436062460229633428422536373014978545251585715765523999626897145245442686109681053275710369838774705588337309829622003283045093978444359731445140597337593964250873138081266774967152048772255402727067611679080157130189967995998175559036833194024776309971929261058416278658972582841062606599798273409641275608450341917573868433271832671690077444780385449073553422161933216384683216622982317479233488903563505597263456144403253067578054817554552829561347839958925083869536682307196235440701185373357913481953704617618446181220385427221091372427837808950835111324692924861139724884562284189627629669549533418733062442369593806837235861863804283010342172306883655302995511580660111634 +/- 5.95e-1485]

```
"""

expected_number_of_cluster_2PD(n::N, θ::T, σ::T, ntrunc::N) where {T<:Number, N<:Integer} = Pkn_2PD_arb.(1:ntrunc, n, θ, σ) |> ar -> map(*, ar, 1:ntrunc) |> sum
expected_number_of_cluster_2PD(n::N, θ::T, σ::T) where {T<:Number, N<:Integer} = expected_number_of_cluster_2PD(n, θ, σ, n)

expected_number_of_cluster_stable(n::N, σ::T, ntrunc::N) where {T<:Number, N<:Integer} = expected_number_of_cluster_2PD(n, 0., σ, ntrunc)
expected_number_of_cluster_stable(n, σ) = expected_number_of_cluster_stable(n, σ, n)


function expected_number_of_clusters_Dirichlet(n::Int64, theta::Float64, ntrunc::Int64)
    return Pkn_Dirichlet_arb.(1:ntrunc, n, theta) |> ar -> map(*, ar, 1:ntrunc) |> sum
end

expected_number_of_clusters_Dirichlet(n::Int64, theta::Float64) = expected_number_of_clusters_Dirichlet(n, theta, n)


function expected_number_of_clusters_PD_direct(n::N, θ::T, σ::T)::T where {T<:Number, N<:Integer}
    return (θ/σ)*(prod([(θ + σ + i -1)/(θ + i -1) for i in 1:n]) -1)
end

function variance_number_of_clusters_PD_direct(n::N, θ::T, σ::T)::T where {T<:Number, N<:Integer}
           if n==1
               return 0
           else
               E_prev =  expected_number_of_clusters_PD_direct(n-1, θ, σ)
               add_term  =  (E_prev*((n-1)*σ - θ*σ) + (n-1)*θ -  σ*σ*(E_prev^2))/ (n- 1 + θ)^2
               return    variance_number_of_clusters_PD_direct(n-1, θ, σ)*(n - 1 + θ + 2*σ)/ (n - 1 + θ) + add_term
           end
       end

function expected_number_of_clusters_Dirichlet_direct(n::N, theta::T) where {T<:Number, N<:Integer}
   return sum( RR(theta)// (RR(theta + i)) for i in 0:(n-1))
end

function variance_number_of_clusters_Dirichlet_direct(n::N, theta::T) where {T<:Number, N<:Integer}
   return sum((RR(theta)*RR(i))//((RR(theta +i))^2) for i in 0:(n-1))
end

function expected_number_of_clusters_Dirichlet_Multinomial_direct(n::N, H::N, theta::T) where {T<:Number, N<:Integer}
    return RR(H) - ((RR(H) -1)*exp( log(risingfac(RR(theta + 1) - RR(theta)//RR(H),n-1)) - log(risingfac(RR(theta + 1),n-1))))
end

expected_number_of_clusters_Dirichlet_Multinomial(n::Int64, H::Int64, theta::Float64) = expected_number_of_clusters_Dirichlet_Multinomial_direct(n, H, theta)


function variance_number_of_clusters(n, exp_nk, pk)
    x = ((1:n).-exp_nk).^2
    v_nk= pk |> ar -> map(*, ar, x) |> sum
    return v_nk
end

variance_number_of_clusters_Dirichlet(n::N, theta::T) where {T<:Number, N<:Integer} = variance_number_of_clusters(n, expected_number_of_clusters_Dirichlet(n, theta), Pkn_Dirichlet_arb.(1:n, n, theta))

variance_number_of_clusters_PD(n::N, theta::T, σ::T) where {T<:Number, N<:Integer} = variance_number_of_clusters(n, expected_number_of_cluster_2PD(n, theta, σ, n) , Pkn_2PD_arb.(1:n, n, theta, σ))


#aliases


"""
    E_2PD(n::N, θ::T, σ::T, ntrunc::N) where {T<:Number, N<:Integer}

Compute the expected number of clusters for a 2-parameter Poisson Dirichlet (also known as Pitman-Yor) prior, by direct computation up to ntrunc (optional argument). If ntrunc = n, this is an exact computation, if it turns out to be too long and that you expect large values of k to have a very small mass, you can set ntrunc < n.

See also: [`E_stable`](@ref), [`E_Dirichlet`](@ref)


# Examples
```julia-repl
julia> E_2PD(10, 0., 0.4)
[2.797180108800000146227991769541811614624157201057103170680445069708107636422615452406497920159297853634382963614055479173158069186418337497999385184360757154429265052398107259805621680470326929725738226729222577234589857779643629658552208848284173508026015152904404059960898279167641852826517393043262627617296658490156505882278711266617066316092443733068941404176241441191328628518109412653612525798394207094225161730723863943059051559749653061842655787927469646092504262924194335937500000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000 +/- 5.60e-1504]

```
"""
E_2PD = expected_number_of_cluster_2PD
"""
    E_PY(n::N, θ::T, σ::T, ntrunc::N) where {T<:Number, N<:Integer}

Compute the expected number of clusters for a 2-parameter Poisson Dirichlet (also known as Pitman-Yor) prior, by direct computation up to ntrunc (optional argument). If ntrunc = n, this is an exact computation, if it turns out to be too long and that you expect large values of k to have a very small mass, you can set ntrunc < n.

See also: [`E_stable`](@ref), [`E_Dirichlet`](@ref)


# Examples
```julia-repl
julia> E_PY(10, 0., 0.4)
[2.797180108800000146227991769541811614624157201057103170680445069708107636422615452406497920159297853634382963614055479173158069186418337497999385184360757154429265052398107259805621680470326929725738226729222577234589857779643629658552208848284173508026015152904404059960898279167641852826517393043262627617296658490156505882278711266617066316092443733068941404176241441191328628518109412653612525798394207094225161730723863943059051559749653061842655787927469646092504262924194335937500000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000 +/- 5.60e-1504]
```
"""
E_PY = expected_number_of_cluster_2PD

"""
    E_stable(n::N, σ::T, ntrunc::N) where {T<:Number, N<:Integer}

Compute the expected number of clusters for a stable process prior, by direct computation up to ntrunc (optional argument). If ntrunc = n, this is an exact computation, if it turns out to be too long and that you expect large values of k to have a very small mass, you can set ntrunc < n.

See also: [`E_2PD`](@ref), [`E_PY`](@ref), [`E_Dirichlet`](@ref)


# Examples
```julia-repl
julia> E_stable(10, 0.1)
[1.317283550395312519640862117090075082018602293487207173467704662870901710848586072338941448910368160445217137507034833714840735500227050610150306349069777704715852016445662205119280375270463908555026018111545594765430207866848367807193962524362807035164072216472928669667852836537566353834476540485547750703723020476001586935706322137304949274062823739150683831653442615917431013917124081952539061271081947740823646018757470139485514508444060324299094553588058194998344774262477585580199956893920898437500000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000 +/- 3.11e-1504]
```
"""
E_stable = expected_number_of_cluster_stable


"""
    E_Dirichlet(n::Int64, theta::Float64, ntrunc::Int64) where {T<:Number, N<:Integer}

Compute the expected number of clusters for a Dirichlet process prior, by direct computation up to ntrunc (optional argument). If ntrunc = n, this is an exact computation, if it turns out to be too long and that you expect large values of k to have a very small mass, you can set ntrunc < n.

See also: [`E_2PD`](@ref), [`E_PY`](@ref), [`E_Dirichlet`](@ref)


# Examples
```julia-repl
julia> E_Dirichlet(10, 0.1)
[1.268596929750522509742103226956398537873085771068571987249592902594795125483947571993880681026641323980721322458062759352392652070330423992468334489803280561487543428122167951021050245266305309523058564612943664653252207376397838984299562945184833920280935894433175454443126982970639290683556142832194488231242118732857777980590412149473759072393894727632777788161259891977237151762150827584780251361172323905738003873059388785341902756850617870450223561134319542176310553349557774240878171574220754051744934444300894535698176817305606825219444113035080791544948019245009777904858615486457170596230701700442717181550577165517662854019961542802495771529879961515496294626049466523909595586762943163930430217545909261330315605742950990173838708604216021274909891705176168468804148449172143254354356340463428016703770246922776019535931638474836624991339643912739739155983869915380851777840945184932873827598663416040011548041974534069037590952150058854019536415499847894848546456328154482639606353077562701342572112539079860044162733170338151051198333830112116544872721287500485113210068274289193032913781410170370106620574401702474340286414131137704367213615618055699430142941143725875436243075702238831256061567411589207476227140330551796033541896922138410130961176340540611096374126740677494288946208178976360070364810387616678096999963421007697822481144326548283583125899364975944100495982998923214048037871138961099375953420318323202550711433054734568046848016758259884095392920728021289115618089199446 +/- 3.29e-1504]
```
"""
E_Dirichlet = expected_number_of_clusters_Dirichlet

"""
E_Dirichlet_multinomial(n::Int64, H::Int64, theta::Float64)

Compute the expected number of clusters for a Dirichlet multinomial process prior, by direct computation.

See also: [`E_2PD`](@ref), [`E_PY`](@ref), [`E_Dirichlet`](@ref)


# Examples
```julia-repl
julia> E_Dirichlet_multinomial(10, 100, 0.1)
[1.2656198162956772044281139501986927560433432728249122838425730172917885456449577827427008029228648259324405450086120700424973872000807961506603125307279025537242800655018045456602670686498182692474826633664464115235080215898472900952036550627041716785370528159699434712650478858211550875546889849261932677562210539024382120856358213274564755169450936766976570333960375938859697028794265405113640006664727462826608028138839169641099450481107571846157530356973423234055102728639482073500720704136233393619357481270310535710377065205775169940077201221510804987245517636551404267972023813211580038972156633810136183697393064520817775221063908711964263522664950624308457120306025925387379275923942136270897436033458705818234588085139913056512529232468622444652256833999355158680011588243488293212125219256186831657014880823083507939645291608505348675502280911274616411712947835536355953711660880159199373720426086123356305436213390603071455111185650982260902472122664033660119541657216813260172107093430967366682626224504543299165643378000241500081885141581110118882802768709028242369726903819960113281002024998943691846080945815486592637086008328168328728610378061578868271472321782307974248638152860658274353781599750807686673500116117928362282072347758716494140304291667565972106364885408051864566662844174074858554185933039304044771000481636810876029924186732368911110789843512235516072634447141633207673761126531255242098046844231863703453219915420808518114669566034153258871751844111589780197562620084 +/- 3.81e-1502]
```
"""
E_Dirichlet_multinomial = expected_number_of_clusters_Dirichlet_Multinomial


"""
   V_PY(n::Int64, theta::Float64, sigma::Float64)
Compute the variance for the number of clusters for σ, by recursion.
See also: [`E_2PD`](@ref), [`E_PY`](@ref), [`E_Dirichlet`](@ref)
# Examples
```julia-repl
julia> V_PY(10, 1.0, 0.2)
2.2826090179905956
```
"""

V_PY  = variance_number_of_clusters_PD_direct


"""
   V_Dirichlet(n::Int64, theta::Float64)
Compute the variance for the number of clusters for Dirichlet process prior.
See also: [`E_2PD`](@ref), [`E_PY`](@ref), [`E_Dirichlet`](@ref)
# Examples
```julia-repl
julia> V_Dirichlet(10, 1.0)
[1.379200522801713277903754094230284706475182665658856135046611237087427563618039808515998992189468379944570420760896951373141849332325522801713277903754094230284706475182665658856135046611237087427563618039808515998992189468379944570420760896951373141849332325522801713277903754094230284706475182665658856135046611237087427563618039808515998992189468379944570420760896951373141849332325522801713277903754094230284706475182665658856135046611237087427563618039808515998992189468379944570420760896951373141849332325522801713277903754094230284706475182665658856135046611237087427563618039808515998992189468379944570420760896951373141849332325522801713277903754094230284706475182665658856135046611237087427563618039808515998992189468379944570420760896951373141849332325522801713277903754094230284706475182665658856135046611237087427563618039808515998992189468379944570420760896951373141849332325522801713277903754094230284706475182665658856135046611237087427563618039808515998992189468379944570420760896951373141849332325522801713277903754094230284706475182665658856135046611237087427563618039808515998992189468379944570420760896951373141849332325522801713277903754094230284706475182665658856135046611237087427563618039808515998992189468379944570420760896951373141849332325522801713277903754094230284706475182665658856135046611237087427563618039808515998992189468379944570420760896951373141849332325522801713277903754094230284706475182665658856135046611237087427563618039808515998992189468379944570420760896951 +/- 4.01e-1504]
```
"""

V_Dirichlet  = variance_number_of_clusters_Dirichlet_direct
