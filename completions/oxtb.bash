# bash completion for oxtb(1).
# generated from src/xhelp.f90 by completions/tools/gen_completions.py, do not edit.
#
# Source of truth: src/xhelp.f90.  Install per user:
#   source ~/.local/share/bash-completion/completions/oxtb
# The file is self contained - it uses only the `compgen` and `compopt`
# builtins, so the bash-completion package (and its _init_completion) is
# not required.  Only the space separated value syntax is offered:
# `--gfn 1` is parsed, `--gfn=1` is silently ignored by the program.

_oxtb_spellings='
  --scc
  --sp
  --grad
  -o
  --opt
  --cycles
  --hess
  --ohess
  --bhess
  --o1nh
  --imagmin
  --imagmax
  --md
  --metadyn
  --bias-input
  --gesc
  --omd
  --metaopt
  --path
  --modef
  --screen
  --mdopt
  --ewin
  --ljson
  --lopt
  --lscan
  --lts
  --lirc
  --lvib
  --method
  --lopt-engine
  --lopt-level
  --lopt-hessian
  --lopt-max-steps
  --lopt-recalc
  --lopt-hess-step
  --lopt-sccacc
  --lopt-hess-jobs
  --lopt-restart
  --lopt-hessian-file
  --lscan-coord
  --lscan-atoms
  --lscan-step
  --lscan-nstep
  --lscan-coord2
  --lscan-atoms2
  --lscan-step2
  --lscan-nstep2
  --lscan-mode
  --lscan-method
  --lscan-level
  --lscan-ndim
  --lscan-restart
  --lscan-oneway
  --lscan-cold-start
  --lscan-fc
  --lscan-expo
  --lscan-max-steps
  --lscan-hess-step
  --lscan-sccacc
  --lscan-hess-jobs
  --lscan-hessian
  --lscan-hessian-file
  --lts-method
  --lts-trust
  --lts-trust-max
  --lts-hlow
  --lts-max-steps
  --lts-level
  --lts-mode-overlap-min
  --lts-alpha
  --lts-micro-cycles
  --lts-target-mode
  --lts-init-mode
  --lts-escape-step
  --lts-delta
  --lts-dr
  --lts-rot-max-iter
  --lts-rot-alpha
  --lts-rot-frms
  --lts-ninit
  --lts-seed
  --lts-step-max
  --lts-step0
  --lts-translate-max-iter
  --lts-recalc
  --lts-no-verify
  --lts-no-rigid
  --lts-restart
  --lts-hessian
  --lts-hessian-file
  --lts-hess-step
  --lts-step
  --lts-sccacc
  --lts-hess-jobs
  --lts-product
  --lts-images
  --lts-kmin
  --lts-kmax
  --lts-kdecay
  --lts-no-dynamic-k
  --lts-climb
  --lts-climb-switch
  --lts-no-idpp
  --lts-no-align
  --lts-relax-ends
  --lts-max-displ
  --lts-frms-avg
  --lts-frms-max
  --lts-autoneb-nstart
  --lts-autoneb-max-paths
  --lts-autoneb-max-depth
  --lts-autoneb-restart
  --lts-string-nnodes
  --lts-string-grow-tol
  --lts-string-ds
  --lts-string-nu2
  --lts-string-grow-steps
  --lts-no-string-spline
  --lts-string-history
  --lirc-method
  --lirc-step
  --lirc-max-steps
  --lirc-level
  --lirc-swap
  --lirc-oneway
  --lirc-target-mode
  --lirc-micro-cycles
  --lirc-micro-thresh
  --lirc-euler-n
  --lirc-loose-cycles
  --lirc-mbs-tol
  --lirc-mbs-max-k
  --lirc-mbs-points
  --lirc-dwi-n
  --lirc-e-rise
  --lirc-erise-steps
  --lirc-grad-tol
  --lirc-osc-window
  --lirc-adaptive
  --lirc-recalc
  --lirc-no-endpoint-relax
  --lirc-endpoint-relax
  --lirc-rmsd-end
  --lirc-hessian
  --lirc-hessian-file
  --lirc-hess-step
  --lirc-step-hess
  --lirc-sccacc
  --lirc-hess-jobs
  --lirc-restart
  --lvib-hessian
  --lvib-hessian-file
  --lvib-hess-step
  --lvib-step
  --lvib-sccacc
  --lvib-hess-jobs
  --lvib-restart
  --lvib-odlr
  --lvib-no-project
  --lvib-temp
  --lvib-rotor-cutoff
  --lvib-imag-cutoff
  --lvib-freq-scal
  --ceh
  --gfn
  --gfn0
  --gfn1
  --gfn2
  --gfnff
  --gff
  --mcgfnff
  --iff
  --nat
  --tblite
  --dipro
  --ptb
  --spinpol
  --vparam
  -c
  --chrg
  --charge
  -u
  --uhf
  --efield
  --etemp
  -a
  --acc
  --iterations
  --maxiterations
  --tm
  --orca
  --mopac
  --driver
  --oniom
  --cut
  --alpb
  -g
  --gbsa
  --gbe
  --gb
  --cosmo
  --tmcosmo
  --cpcmx
  --cma
  --esp
  --stm
  --pop
  --molden
  --dipole
  --raman
  --alpha
  --wbo
  --lmo
  --fod
  --enso
  --vip
  --vea
  --vipea
  --vfukui
  --vomega
  --color
  --namespace
  --json
  --help-values
  --wrtopo
  --citation
  --license
  --version
  -v
  --verbose
  -V
  --very-verbose
  --ceasefiles
  --strict
  -h
  --help
  --coffee
  -I
  --input
  --copy
  --nocopy
  --restart
  --norestart
  --define
  --nocellopt
  -P
  --parallel
  --blas-threads
'

# file names: the geometry argument and the file taking options
_oxtb_files()
{
    COMPREPLY=( $(compgen -f -- "$1") )
    if [ ${#COMPREPLY[@]} -gt 0 ] && type compopt >/dev/null 2>&1; then
        compopt -o filenames 2>/dev/null
    fi
}

_oxtb()
{
    local cur prev
    cur=${COMP_WORDS[COMP_CWORD]}
    prev=''
    if [ "${COMP_CWORD}" -gt 0 ]; then
        prev=${COMP_WORDS[COMP_CWORD-1]}
    fi
    case $prev in
    -o|--opt)
        COMPREPLY=( $(compgen -W 'crude sloppy loose lax normal tight verytight vtight extreme superloose' -- "$cur") )
        return 0
        ;;
    --ohess)
        COMPREPLY=( $(compgen -W 'crude sloppy loose lax normal tight verytight vtight extreme superloose' -- "$cur") )
        return 0
        ;;
    --bhess)
        COMPREPLY=( $(compgen -W 'crude sloppy loose lax normal tight verytight vtight extreme superloose' -- "$cur") )
        return 0
        ;;
    --metaopt)
        COMPREPLY=( $(compgen -W 'crude sloppy loose lax normal tight verytight vtight extreme superloose' -- "$cur") )
        return 0
        ;;
    --method)
        COMPREPLY=( $(compgen -W 'auto sd cg lbfgs bfgs diis rfo prfo bprfo dimer neb cineb autoneb string gsm gs hpc lqa eulerpc ancopt rf lancopt l-ancopt fire inertial pbc_lbfgs' -- "$cur") )
        return 0
        ;;
    --lopt-engine)
        COMPREPLY=( $(compgen -W 'sd cg lbfgs bfgs diis rfo' -- "$cur") )
        return 0
        ;;
    --lopt-level)
        COMPREPLY=( $(compgen -W 'superloose loose normal tight vtight' -- "$cur") )
        return 0
        ;;
    --lopt-hessian)
        COMPREPLY=( $(compgen -W 'numerical model file' -- "$cur") )
        return 0
        ;;
    --lopt-restart)
        COMPREPLY=( $(compgen -W 'off pin warm' -- "$cur") )
        return 0
        ;;
    --lscan-coord)
        COMPREPLY=( $(compgen -W 'distance angle dihedral' -- "$cur") )
        return 0
        ;;
    --lscan-mode)
        COMPREPLY=( $(compgen -W 'relaxed rigid' -- "$cur") )
        return 0
        ;;
    --lscan-method)
        COMPREPLY=( $(compgen -W 'sd cg lbfgs bfgs diis rfo' -- "$cur") )
        return 0
        ;;
    --lscan-level)
        COMPREPLY=( $(compgen -W 'superloose loose normal tight vtight' -- "$cur") )
        return 0
        ;;
    --lscan-ndim)
        COMPREPLY=( $(compgen -W '1 2' -- "$cur") )
        return 0
        ;;
    --lscan-restart)
        COMPREPLY=( $(compgen -W 'off pin warm' -- "$cur") )
        return 0
        ;;
    --lscan-hessian)
        COMPREPLY=( $(compgen -W 'numerical model file' -- "$cur") )
        return 0
        ;;
    --lts-method)
        COMPREPLY=( $(compgen -W 'prfo bprfo dimer neb cineb autoneb string gsm' -- "$cur") )
        return 0
        ;;
    --lts-level)
        COMPREPLY=( $(compgen -W 'superloose loose normal tight vtight' -- "$cur") )
        return 0
        ;;
    --lts-init-mode)
        COMPREPLY=( $(compgen -W 'h0 force' -- "$cur") )
        return 0
        ;;
    --lts-ninit)
        COMPREPLY=( $(compgen -W 'random force' -- "$cur") )
        return 0
        ;;
    --lts-restart)
        COMPREPLY=( $(compgen -W 'off pin warm' -- "$cur") )
        return 0
        ;;
    --lts-hessian)
        COMPREPLY=( $(compgen -W 'numerical model file' -- "$cur") )
        return 0
        ;;
    --lirc-method)
        COMPREPLY=( $(compgen -W 'gs hpc lqa eulerpc' -- "$cur") )
        return 0
        ;;
    --lirc-level)
        COMPREPLY=( $(compgen -W 'superloose loose normal tight vtight' -- "$cur") )
        return 0
        ;;
    --lirc-hessian)
        COMPREPLY=( $(compgen -W 'numerical model file' -- "$cur") )
        return 0
        ;;
    --lirc-restart)
        COMPREPLY=( $(compgen -W 'off pin warm' -- "$cur") )
        return 0
        ;;
    --lvib-hessian)
        COMPREPLY=( $(compgen -W 'numerical model file' -- "$cur") )
        return 0
        ;;
    --lvib-restart)
        COMPREPLY=( $(compgen -W 'off pin warm' -- "$cur") )
        return 0
        ;;
    --gfn)
        COMPREPLY=( $(compgen -W '0 1 2 ff' -- "$cur") )
        return 0
        ;;
    --alpb)
        COMPREPLY=( $(compgen -W 'acetone acetonitrile aniline benzaldehyde benzene ch2cl2 chcl3 cs2 dioxane dmf dmso ethanol ether ethylacetate furane hexadecane hexane methanol nitromethane octanol woctanol phenol toluene thf water' -- "$cur") )
        return 0
        ;;
    -g|--gbsa)
        COMPREPLY=( $(compgen -W 'acetone acetonitrile benzene ch2cl2 chcl3 cs2 dmf dmso ether water methanol hexane thf toluene' -- "$cur") )
        return 0
        ;;
    --color)
        COMPREPLY=( $(compgen -W 'auto always never' -- "$cur") )
        return 0
        ;;
    --wrtopo)
        COMPREPLY=( $(compgen -W 'etot gnorm nb bpair alist blist tlist vtors vbond vangl hbbond eeq json' -- "$cur") )
        return 0
        ;;
    --bias-input|--gesc)
        _oxtb_files "$cur"
        return 0
        ;;
    --path)
        _oxtb_files "$cur"
        return 0
        ;;
    --lopt-hessian-file)
        _oxtb_files "$cur"
        return 0
        ;;
    --lscan-hessian-file)
        _oxtb_files "$cur"
        return 0
        ;;
    --lts-hessian-file)
        _oxtb_files "$cur"
        return 0
        ;;
    --lirc-hessian-file)
        _oxtb_files "$cur"
        return 0
        ;;
    --lvib-hessian-file)
        _oxtb_files "$cur"
        return 0
        ;;
    --vparam)
        _oxtb_files "$cur"
        return 0
        ;;
    --driver)
        _oxtb_files "$cur"
        return 0
        ;;
    -I|--input)
        _oxtb_files "$cur"
        return 0
        ;;
    --cycles|--imagmin|--imagmax|--metadyn|--modef|--ewin|--lopt-max-steps|--lopt-recalc|--lopt-hess-step|--lopt-sccacc|--lopt-hess-jobs|--lscan-atoms|--lscan-step|--lscan-nstep|--lscan-fc|--lscan-expo|--lscan-max-steps|--lscan-hess-step|--lscan-sccacc|--lscan-hess-jobs|--lts-trust|--lts-trust-max|--lts-hlow|--lts-max-steps|--lts-mode-overlap-min|--lts-alpha|--lts-micro-cycles|--lts-target-mode|--lts-escape-step|--lts-delta|--lts-dr|--lts-rot-max-iter|--lts-rot-alpha|--lts-rot-frms|--lts-seed|--lts-step-max|--lts-step0|--lts-translate-max-iter|--lts-recalc|--lts-hess-step|--lts-step|--lts-sccacc|--lts-hess-jobs|--lts-product|--lts-images|--lts-kmin|--lts-kmax|--lts-kdecay|--lts-climb-switch|--lts-max-displ|--lts-frms-avg|--lts-frms-max|--lts-autoneb-nstart|--lts-autoneb-max-paths|--lts-autoneb-max-depth|--lts-autoneb-restart|--lts-string-nnodes|--lts-string-grow-tol|--lts-string-ds|--lts-string-nu2|--lts-string-grow-steps|--lirc-step|--lirc-max-steps|--lirc-target-mode|--lirc-micro-cycles|--lirc-micro-thresh|--lirc-euler-n|--lirc-loose-cycles|--lirc-mbs-tol|--lirc-mbs-max-k|--lirc-mbs-points|--lirc-dwi-n|--lirc-e-rise|--lirc-erise-steps|--lirc-grad-tol|--lirc-osc-window|--lirc-recalc|--lirc-rmsd-end|--lirc-hess-step|--lirc-step-hess|--lirc-sccacc|--lirc-hess-jobs|--lvib-hess-step|--lvib-step|--lvib-sccacc|--lvib-hess-jobs|--lvib-temp|--lvib-rotor-cutoff|--lvib-imag-cutoff|--lvib-freq-scal|--ceh|--nat|--dipro|-c|--chrg|--charge|-u|--uhf|--efield|--etemp|-a|--acc|--iterations|--maxiterations|--oniom|--gbe|--gb|--cosmo|--tmcosmo|--cpcmx|--namespace|-P|--parallel|--blas-threads)
        # an argument without a value set: no candidates
        COMPREPLY=()
        return 0
        ;;
    esac
    if [[ $cur == -* ]]; then
        COMPREPLY=( $(compgen -W "$_oxtb_spellings" -- "$cur") )
        return 0
    fi
    _oxtb_files "$cur"
}

complete -F _oxtb oxtb
