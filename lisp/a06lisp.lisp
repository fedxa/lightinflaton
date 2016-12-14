;; Derived from fortan/a06.f
;;    DPDFS are not implemented!
;; These are AMP06 set from Alekhin's NNLO
;; The AMP06 set -- the NNLO PDFs extracted from the the inclusive DIS and
;;                  fixed-target Drell-Yan data
;;                  [hep-ph/0606237] Phys. Rev. D74, 054033 (2006)
;;    http://mail.ihep.ru/~alekhin/pdfa06/amp06.html
;; (see http://durpdg.dur.ac.uk/hepdata/pdf.html for other PDFs)
;;
;;     This is a code for the NNLO parton distributions 
;;     in the variable-flavor-number (VFN) schem with account 
;;     of their experimental and theoretical uncertainties. 
;;     The Q**2 range is 0.8d0 < Q**2 < 2d8, the x range is 1d-7 < x < 1d0
;;     (for the values of PDFs and strong coupling constant at Q2 < 0.8 GeV^2 
;;     (x < 1d-7) their values at Q^2 = 0.8 GeV^2 (x = 1d-7) are returned).
;;
;;  Output parameters:
;;     The array PDFS contains fitted values of the strong coupling constant 
;;     and the parton distributions at given x and Q:
;;        PDFS(0) -- \alpha_s
;;        PDFS(1) -- valence u-quarks 
;;        PDFS(2) -- valence d-quarks
;;        PDFS(3) -- gluons 
;;        PDFS(4) -- sea u-quarks 
;;        PDFS(5) -- s-quarks 
;;        PDFS(6) -- sea d-quarks 
;;        PDFS(7) -- c-quarks
;;        PDFS(8) -- b-quarks
;;        PDFS(9) -- t-quarks
;;     NPDF is the number of PDFs returned (NPDF=9 for the VFN scheme).
;;     Output array DPDFS(0:npdf,npar) contains derivatives of \alpha_s and
;;     the PDFs on the fitted parameters with the number of the parameters 
;;     returned in NPAR. With the derivatives of \alpha_s included one can take 
;;     into account the correlations of the fitted PDFs with \alpha_s as well.
;;     All derivatives are transformed to the orthonormal 
;;     basis of eigenvectors of the parameters error matrix. For this reason 
;;     the variation of the PDFs in the derivatives directions can be performed 
;;     independently. For example the dispersion of the i-th PDF can be stored 
;;     in DELPDF using the code 

(defpackage :a06
  (:use :common-lisp :gsl)
  (:export "A06-INIT" "A06"))

(in-package :a06)

(eval-when (:compile-toplevel :load-toplevel :execute)
  (unless (fboundp 'concat)
    (defmacro concat (&rest args)
      `(concatenate 'string ,@args))))

(defconstant +nxb+ 99)
(defconstant +nxbb+ (truncate +nxb+ 2))
(defconstant +nq+ 20)
(defconstant +np+ 9)
(defconstant +nvar+ 23)
(defconstant +npar+ +nvar+)

(defparameter *locdir* "./")
(defparameter +nexp+ #(0 3 4 5 5 5 5 5 5 5))
(defconstant +xmin+ 1d-7)
(defconstant +xmax+ 1d0)
(defconstant +qsqmin+ 0.8d0)
(defconstant +qsqmax+ 2d8)
(defconstant +dels+ (/ (- (log (log (/ +qsqmax+ 0.04d0)))
			  (log (log (/ +qsqmin+ 0.04d0))))
		       (1- +nq+)))

(defconstant +npdf+ 9)

(defconstant +x1+ 0.3d0)
(defconstant +xlog1+ (log +x1+))
(defconstant +delx+ (/ (- (log +x1+) (log +xmin+))
		       (1- +nxbb+)))
(defconstant +DELX1+ (/ (expt (- 1.d0 +x1+) 2)
			(1+ +nxbb+) ))

(defvar *f*)
(defvar *fsplines*)
(defvar *df*)

;; X GRID
(defvar *xx*)
(defmacro *xx* (i) `(grid:aref *xx* ,i))
(defun xx-grid-init ()
  (setq *xx* (grid:make-foreign-array 'double-float :dimensions +nxb+))
  (loop :for kx :from 0 :below +nxbb+ :do
     (setf (*xx* kx) (exp (+ (log +xmin+)
				  (* +delx+ kx)))))
  (loop :for kx :from +nxbb+ :below (1- +nxb+) :do
     (setf (*xx* kx) (- 1 (sqrt (abs (- (expt (- 1 +x1+) 2)
					     (* +delx1+ (- kx +nxbb+ -1))))))))
  (setf (*xx* (1- +nxb+)) 1.0d0))

(defmacro *f* (n m i) `(grid:aref (aref *f* ,i ,m) ,n))
(defun read-pdfs ()
  ;; Read input tables
  (print "***** Reading PDFs from tables *****")
  (setq *f* (make-array (list (1+ +npdf+) +nq+)))
  (loop :for i :from 0 :to +npdf+ :do
     (loop :for m :from 0 :below +nq+ :do
	(setf (aref *f* i m)
	      (grid:make-foreign-array 'double-float :dimensions +nxb+))))
  (with-open-file (nport (concat *locdir* "a06.pdfs_3_vfn"))
    (loop :for n :from 0 :below (1- +nxb+) :do
       (loop :for m :from 0 :below +nq+ :do
	  (loop :for i :from 0 :to +npdf+ :do
	     (setf (grid:aref (aref *f* i m) n)
		   (coerce (read nport) 'double-float)))))
    (loop :for i :from 0 :to +npdf+ :do
       (loop :for m :from 0 :below +nq+ :do
	  (if (/= i 0)
	      (setf (*f* (1- +nxb+) m i ) 0d0)
	      (setf (*f* (1- +nxb+) m i ) (*f* (- +nxb+ 2) m i)))
	  (loop :for n :from 0 :below (1- +nxb+) :do
	     (setf (*f* n m i) (/ (*f* n m i)
				  (expt (- 1 (*xx* n)) (elt +nexp+ i)))))))))

(defun a06-init ()
  (xx-grid-init)
  (read-pdfs)
  (setq *fsplines* (make-array (list (1+ +npdf+) +nq+)))
  (loop :for i :from 0 :to +npdf+ :do
     (loop :for m :below +nq+ :do
	(setf (aref *fsplines* i m)
	      (make-spline +cubic-spline-interpolation+ *xx* (aref *f* i m))))))

(a06-init)

(defun a06 (xb q2 i)
  (if (< q2 +qsqmin+) (progn
			(setq q2 +qsqmin+)
			(print "q2 below minimum")))
  (if (> q2 (* 0.99 +qsqmax+)) (progn
			(setq q2 (* 0.99 +qsqmax+))
			(print "q2 above maximum")))
  (if (< xb +xmin+) (progn
			(setq xb +xmin+)
			(print "xb below minimum")))
  (if (> xb +xmax+) (progn
			(setq xb +xmax+)
			(print "xb above maximum")))
  (let* ((ss (- (log (log (/ q2 0.04d0))) (log (log (/ +qsqmin+ 0.04d0)))))
	 (m (truncate ss +dels+))
	 (b (- (/ ss +dels+) m))
	 f0 fp fm pdf)
    (setq f0 (evaluate (aref *fsplines* i m) xb))
    (setq fp (evaluate (aref *fsplines* i (1+ m)) xb))
    (if (>= m 1)
	(progn
	  (setq fm (evaluate (aref *fsplines* i (1- m)) xb))
	  (setq pdf (+ (/ (* fm b (- b 1)) 2)
		       (* f0 (- 1 (* b b)))
		       (/ (* fp b (+ b 1)) 2))))
	(setq pdf (+ (* f0 (- 1 b)) (* fp b))))
    (* pdf (expt (- 1 xb) (elt +nexp+ i)))))

;; Seems to be even better, than the original
