;; Should be actually moved into maxima code for convenience...
;; Not much of a speed critical things here

(eval-when (:compile-toplevel :load-toplevel :execute)
  (ql:quickload "gsll")
  (load "a06lisp")
  )

(defpackage :chidecays
  (:use :common-lisp :gsll :a06))

(in-package :chidecays)

(eval-when (:compile-toplevel :load-toplevel :execute)
  (unless (fboundp 'concat)
    (defmacro concat (&rest args)
      `(concatenate 'string ,@args))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
(defmacro defmnumfun (fun)
  "This redefines the function as a maxima totally numerical function.
If maxima is not loaded it creates a dummy :maxima package."
  (if (not (find-package :maxima))
      (make-package :maxima))
  (let* ((f fun)
	 (mfun (intern (concat "$" (symbol-name f)) (find-package :maxima)))
	 (args (gensym)))
    `(defun ,mfun (&rest ,args)
       (if (every #'numberp ,args)
	   (apply #',f ,args)
	   (cons '(,mfun) ,args)))))
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;All masses are in GeV
(defparameter +GF+ 1.16637d-5)
(defparameter +me+ 0.510998910d-3)
(defparameter +mmu+ 105.6583668d-3)
(defparameter +mtau+ 1776.84d-3)
(defparameter +mpi+ 139.57018d-3)
(defparameter +mpi0+ 134.9766d-3)
(defparameter +mW+ 80.398d0)
(defparameter +mZ+ 91.1876)
(defparameter +md+ 5.04d-3)
(defparameter +mu+ 2.55d-3)
(defparameter +ms+ 105d-3)
(defparameter +mc+ 1.27d0)
(defparameter +mb+ 4.20d0)
(defparameter +mt+ 171.3d0)

(defparameter +mK+ 493.677d-3)
(defparameter +mK0+ 497.614d-3)

(defparameter +mDmeson+ 1865d-3)
(defparameter +mBmeson+ 5279d-3)

(defparameter +alpha+ (/ 137.035999679d0))
(defparameter +alphaSMZ+ 0.1176d0)

(defun a06_gluon (xb q2)
  (a06:a06 xb q2 3))
(defun a06_alphaS (q2)
  (a06:a06 0.5d0 q2 0))

(defvar *muF2*)
(defvar *tau*)
;; My calculation
(defun a12 (tau)
  (let ((f (if (<= tau 1)
	       (expt (asin (sqrt tau)) 2)
	       (let ((tt (sqrt (- 1 (/ 1 tau)))))
		 (- (/ (expt (- (log (/ (+ 1 tt) (- 1 tt))) (* PI #C(0 1))) 2) 4))))))
     (/ (* 2 (+ tau (* f (- tau 1))))
	tau tau)))
(defun a12sum (mH)
  (let ((mH4 (* mH mH 0.25)))
    (expt (* 0.75 (abs (+ (a12 (/ mh4 +mc+ +mc+))
			  (a12 (/ mh4 +mb+ +mb+))
			  (a12 (/ mh4 +mt+ +mt+)))))
	 2)))

(defun dldtfunc (x)
  (/ (* (a06_gluon x *muF2*) (a06_gluon (/ *tau* x) *muF2*)) x))

;; This works (no boundery checks!) but is two times slower, than fortran
(defun sigma_2 (mH s rcoeff)
  (let ((*muF2* (* mH mH rcoeff rcoeff))
	(*tau* (/ (* mH mH) s)))
    (/ (* (integration-QAG #'dldtfunc *tau* 1.0 :gauss15
			   -1.0 ; no absolute error
			   1d-3 ; relative error
			   )
	  +GF+ (expt (a06_alphaS *muF2*) 2) (a12sum mH))
       (* 288 (sqrt 2.0) PI))))

(defmnumfun sigma_2)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;Quark MSbar masses
(defun cfunc (mu)
  (let ((x (/ (alphaS mu) PI)))
    (if (< mu +mc+)
	(* (expt (* 9/2 x) 4/9)
	   (+ 1 (* x (+ 0.895 (* x (+ 1.371 (* x 1.952)))))))
	(if (< mu +mb+)
	    (* (expt (* 25/6 x) 12/25)
	       (+ 1 (* x (+ 1.014 (* x (+ 1.389 (* x 1.091)))))))
	    (if (< mu +mt+)
		(* (expt (* 23/6 x) 12/23)
		   (+ 1 (* x (+ 1.175 (* x (+ 1.501 (* x 0.1725)))))))
		(* (expt (* 7/2 x) 4/7)
		   (+ 1 (* x (+ 1.398 (* x (+ 1.793 (* x -0.6834))))))))))))
(defun ms (mu)
  (* 104d-3 (/ (cfunc mu) (cfunc +mtau+))))
(defun mc (mu)
  (* +mc+ (/ (cfunc mu) (cfunc +mc+))))
(defun mb (mu)
  (* +mb+ (/ (cfunc mu) (cfunc +mb+))))
(defmnumfun ms)
(defmnumfun mc)

;Lepton decays
(defun gammall (mchi beta ml)
  (if (> mchi (* 2 ml))
      (/ (* beta ml ml
	    (expt (- 1 (/ (* 4 ml ml) mchi mchi)) 3/2))
	 4 PI mchi)
      0))

(defun gammaee (mchi beta)
  (gammall mchi beta +me+))
(defun gammamumu (mchi beta)
  (gammall mchi beta +mmu+))
(defun gammatautau (mchi beta)
  (gammall mchi beta +mtau+))
(defun gammass (mchi beta)
  (let ((a (/ (alphaS mchi) PI)))
    (* (* 3 (gammall mchi beta (ms mchi)))
       (+ 1 (* a (+ 5.67))))))
(defun gammassKcut (mchi beta)
  (if (> mchi (* 2 +mK+))
      (gammass mchi beta)
      0))
(defun gammass0 (mchi beta)
  (* (* 3 (gammall mchi beta +ms+))))
(defun gammacc (mchi beta)
  (if (> mchi (* 2 +mDmeson+))
      (let ((a (/ (alphaS mchi) PI)))
	(* (* 3 (gammall mchi beta (mc mchi)))
	   (+ 1 (* a (+ 5.67)))))
      0))
(defun gammacc0 (mchi beta)
  (if (> mchi (* 2 +mDmeson+))
      (* 3 (gammall mchi beta +mc+))
      0))
(defun gammabb (mchi beta)
  (if (> mchi (* 2 +mBmeson+))
      (let ((a (/ (alphaS mchi) PI)))
	(* (* 3 (gammall mchi beta (mb mchi)))
	   (+ 1 (* a (+ 5.67)))))
      0))
(defun gammatt (mchi beta)
  (* 3 (gammall mchi beta +mt+)))

;Meson decays
(defun gammaf0f0 (mchi beta mf)
  (if (> mchi (* 2 mf))
      (/ (* beta mchi
	    (expt (+ 2/9 (/ (* 11/9 mf mf) mchi mchi)) 2)
	    (sqrt (- 1 (/ (* 4 mf mf) mchi mchi))))
	 16 pi)
      0))

(defun gammapipi (mchi beta)
  (* 2 (gammaf0f0 mchi beta +mpi+)))
(defun gammapi0pi0 (mchi beta)
  (gammaf0f0 mchi beta +mpi0+))
(defun gammaKKchiral (mchi beta)
  (* 2 (gammaf0f0 mchi beta +mK+)))
(defun gammaK0K0chiral (mchi beta)
  (* 2 (gammaf0f0 mchi beta +mK0+)))


;; These are just for show -- not real formulas!!!!
(defun gammaKK (mchi beta)
  (let ((x (gammaKKchiral mchi beta))
	(y (/ (gammass mchi beta) 2)))
    (/ (* x y) (+ x y))))
(defun gammaK0K0 (mchi beta)
  (let ((x (gammaK0K0chiral mchi beta))
	(y (/ (gammass mchi beta) 2)))
    (/ (* x y) (+ x y))))


(defun make-spline-pairlist (list type)
  (let ((x (grid:make-foreign-array 'double-float :initial-contents
			(loop for f in list collect (first f))))
	(y (grid:make-foreign-array 'double-float :initial-contents
			(loop for f in list collect (second f)))))
    (make-spline type x y)))
(defmacro definterpolation (fname data)
  (let ((func (intern fname))
	(m (gensym))
	(acc (gensym)))
    `(eval-when (:load-toplevel :execute)
       (let* ((,acc (make-acceleration))
	      (spline (make-spline-pairlist ,data +linear-interpolation+)))
	 (defun ,func (,m)
	   (evaluate spline ,m :acceleration ,acc))))))
  

;;Enhanced to pipi decay -- phenomenological
(definterpolation "RATIOPIPIMUMU"
    '(
;      (3.00895E-1	4.58601E+0)
      (3.00895E-1	1.58601E+0)
      (4.19882E-1	7.60993E+0)
      (6.00949E-1	1.21441E+1)
      (7.97556E-1	1.97691E+1)
      (9.42801E-1	7.78151E+1)
      (9.74610E-1	1.85547E+2)
      (9.87942E-1	2.41353E+2)
      (9.95567E-1	2.22743E+2)
;    (9.90010E-1	1.69258E+2)
;    (9.92345E-1	1.34372E+2)
      (1.01006E+0	8.00972E+1)
      (1.05375E+0	4.20843E+1)
      (1.09500E+0	2.42281E+1)
      (1.16214E+0	9.45603E+0)
      (1.24225E+0	1.65230E+0)
      ))
;;Enhanced to KK decay -- phenomenological
(definterpolation "RATIOKKMUMU"
    '((9.86220E-1	1.81770E+0)
      (9.88912E-1	1.65447E+1)
      (9.91682E-1	4.21244E+1)
      (9.96949E-1	5.52994E+1)
      (1.00738E+0	6.76958E+1)
      (1.05910E+0	6.68872E+1)
      (1.09527E+0	6.14375E+1)
      (1.16755E+0	4.35612E+1)
      (1.24247E+0	3.11097E+1)
      (1.32257E+0	2.17555E+1)
      (1.39754E+0	1.70559E+1)
      ;; These are my fantasy -- smooth interpoaltion to ss
      (1.5E+0   	1E+1)
      (2E+0    	0.41E+1)     
      ))

(defun gammapipitot (mchi beta)
  (cond
    ((< mchi 3.00895E-1) (+ (gammapipi mchi beta) (gammapi0pi0 mchi beta)))
    ((< mchi 1.24225E+0) (* (ratiopipimumu (coerce mchi 'double-float)) (gammamumu mchi beta)))
    (t -1)))
;    (t (gammagg mchi beta))))


;; Do not believe this above mchi>1.2!!!!
(defun gammaKKtot (mchi beta)
  (if (<= mchi (* 2 +mK+))
      0
      (if (and (> mchi 9.86220E-1) (< mchi 2))
	  (* (ratioKKmumu (coerce mchi 'double-float)) (gammamumu mchi beta))
	  (+ (gammaKK mchi beta) (gammaK0K0 mchi beta))
	  )))


;Photon decay
(defun xfunc (y)
  (if (> y 1)
      (atan (/ (sqrt (- y 1))))
      (* 1/2 (+ pi (* #C(0 1) (log (/ (+ 1 (sqrt (- 1 y)))
				      (- 1 (sqrt (- 1 y))))))))))
(defun x2func (y)
  (expt (xfunc y) 2))
(defun Fw (y)
  (+ 2 (* 3 y (+ 1 (* (- 2 y) (x2func y))))))
(defun Ff (y)
  (* -2 y (+ 1 (* (- 1 y) (x2func y)))))
(defun Fgamma (mchi)
  (let ((y (lambda (m) (/ (* 4 m m) mchi mchi))))
    (+ (Fw (funcall y +mW+))
       (loop
	  for f in `((1 1 ,+me+)
		     (1 1 ,+mmu+)
		     (1 1 ,+mtau+)
		     (3 -1/3 ,+md+)
		     (3 2/3 ,+mu+)
		     (3 -1/3 ,+ms+)
		     (3 2/3 ,+mc+)
		     (3 -1/3 ,+mb+)
		     (3 2/3 ,+mt+)
		     )
	  sum (* (first f) (expt (second f) 2) (Ff (funcall y (third f))))))))

;; Be cartefull -- it is ok only for small mh! (check at okun') ?
(defun gammagammagamma (mchi beta)
  (let ((F (Fgamma mchi)))
    (/ (* beta mchi +alpha+ +alpha+
	  (realpart (* F (conjugate F))))
       128 pi pi pi)))

(defun alphaSlambda (mu l nf)
  (let* ((b0 (- 11 (* 2/3 nf)))
	 (b1 (- 51 (* 19/3 nf)))
	 (b2 (+ 2857 (* -5033/9 nf) (* 325/27 nf nf)))
	 (lm (* 2 (log (/ mu l))))
	 (llm (log lm)))
    (/ (* 4 PI (+ 1
		   (/ (* -2 b1 llm) b0 b0 lm)
		   (/ (* 4 b1 b1 (+ (expt (- llm 1/2) 2)
				    (/ (* b2 b0) 8 b1 b1)
				    -5/4))
		      (expt b0 4) lm lm)))
       b0 lm)))

;(asdf:oos 'asdf:load-op 'gsll)
;(gsll:make-one-dimensional-root-solver-f gsll:+brent-fsolver+ 'quadratic 0.0d0 5.0d0)

    
(defun alphaSrun (mu mu0 as0 nf)
  (/ as0
     (+ 1 (* as0 (/ (- 33 (* 2 nf)) 6 PI) (log (/ mu mu0))))))
(defparameter +alphaSmt+ (alphaSrun +mt+ +MZ+ +alphaSMZ+ 5))
(defparameter +alphaSmb+ (alphaSrun +mb+ +MZ+ +alphaSMZ+ 5))
(defparameter +alphaSmc+ (alphaSrun +mc+ +mb+ +alphaSmb+ 4))
;; Rough approximation for alphaS at scale mu
;; IMPORTANT: mu is the scale in GeV, while in alphas_f90 the argument is _sqared_ scale!
(defun alphaS (mu)
  (if (< mu +mc+)
      (alphaSrun mu +mc+ +alphaSmc+ 3)
      (if (< mu +mb+)
	  (alphaSrun mu +mb+ +alphaSmb+ 4)
	  (if (< mu +mt+)
	      (alphaSrun mu +mt+ +alphaSmt+ 5)
	      (alphaSrun mu +mt+ +alphaSmt+ 6)))))

(defun Fg (mchi)
  (let ((y (lambda (m) (/ (* 4 m m) mchi mchi))))
    (loop
     for f in (list +md+ +mu+ +ms+ +mc+ +mb+ +mt+)
     sum (* (Ff (funcall y f))))))
(defun gammagg (mchi beta)
  (let ((F (Fg mchi)))
    (/ (* beta mchi (expt (alphaS mchi) 2)
	  (realpart (* F (conjugate F))))
       64 pi pi pi)))

(defun gammatotle (mchi beta)
  (loop
     for g in '(gammaee gammamumu gammatautau
		gammapipitot
		gammaKKtot
		gammagammagamma
		)
     sum (funcall g mchi beta)))

(defun gammatothe (mchi beta)
  (loop
     for g in '(gammaee gammamumu gammatautau
		gammagg gammassKcut
		gammacc gammabb gammatt
		)
     sum (funcall g mchi beta)))

(defun gammatot (mchi beta)
  (if (< mchi 1.23)
      (gammatotle mchi beta)
      (if (< mchi 1.76)
	  ;; Very wild hand made fit in the intermediate region of masses
	  (* 1.1 (gammaKKtot mchi beta))
	  (gammatothe mchi beta))))

(defmacro defbranching (proc)
  (let ((p proc)) 
    `(eval-when (:compile-toplevel :load-toplevel :execute)
      (defun ,(intern (concat "BR" p)) (mchi)
	(/ (,(intern (concat "GAMMA" p)) mchi 1) (gammatot mchi 1)))
      (defmnumfun ,(intern (concat "BR" p))))))
(defbranching "EE")
(defbranching "MUMU")
(defbranching "TAUTAU")
(defbranching "PIPITOT")
(defbranching "KKTOT")
(defbranching "GAMMAGAMMA")


(defbranching "GG")
(defbranching "SS")
(defbranching "SSKCUT")
(defbranching "CC")
(defbranching "BB")
(defbranching "TT")

(defmnumfun GAMMATOT)
(defmnumfun GAMMATOTLE)
(defmnumfun GAMMATOTHE)

(defmnumfun GAMMAEE)
(defmnumfun GAMMAMUMU)
(defmnumfun GAMMAGAMMAGAMMA)
(defmnumfun GAMMAPIPI)
(defmnumfun GAMMAPI0PI0)
(defmnumfun GAMMAPIPITOT)
(defmnumfun GAMMAKKTOT)
(defmnumfun GAMMAGG)
(defmnumfun GAMMASS)
(defmnumfun GAMMASS0)
(defmnumfun GAMMASSKCUT)
(defmnumfun GAMMAKKCHIRAL)
(defmnumfun GAMMAK0K0CHIRAL)
(defmnumfun GAMMAKK)
(defmnumfun GAMMACC)
(defmnumfun GAMMACC0)
(defmnumfun GAMMAK0K0)

(defmnumfun ALPHAS)
