
#ifndef __has_c2_function_hh
#define __has_c2_function_hh 1

#ifdef _MSC_VER
#define _USE_MATH_DEFINES
#define c2_isnan _isnan
#define c2_isfinite _finite
#else
#define c2_isnan std::isnan
#define c2_isfinite std::isfinite
#endif

#include <cmath>
#include <vector>
#include <utility>
#include <string>
#include <stdexcept>
#include <typeinfo>
#include <sstream>
#include <limits> // fails under gcc-4.3 without this here, was ok in c2_function.cc before

/// \brief the exception class for c2_function operations.
class c2_exception : public std::exception {
public:
    /// \brief construct the exception with an error message
    /// \param msgcode the message
    c2_exception(const char msgcode[]) : info(msgcode) { }
    virtual ~c2_exception() throw() { }
    /** Returns a C-style character string describing the general cause
    *  of the current error.  */
    virtual const char* what() const throw() { return info.c_str(); }
private:
    std::string info;
};

// put these forward references here, and with a bogus typename to make swig happy.
template <typename float_type> class c2_composed_function_p;
template <typename float_type> class c2_sum_p;
template <typename float_type> class c2_diff_p;
template <typename float_type> class c2_product_p;
template <typename float_type> class c2_ratio_p;
template <typename float_type> class c2_piecewise_function_p;
template <typename float_type> class c2_quadratic_p;
template <typename float_type> class c2_ptr;

template <typename float_type> class c2_fblock
{
public:
        /// \brief the abscissa
        float_type x;
        /// \brief the value of the function at \a x
        float_type y;
        /// \brief the derivative at \a x
        float_type yp;
        /// \brief the second derivative at \a x
        float_type ypp;
        /// flag, filled in by c2_function::fill_fblock(), indicating the derivative is NaN of Inf
        bool ypbad;
        /// flag, filled in by c2_function::fill_fblock(), indicating the second derivative is NaN of Inf
        bool yppbad;
};

template <typename float_type=double> class c2_function {
public:
    /// \brief get versioning information for the header file
    /// \return the CVS Id string
        const std::string cvs_header_vers() const { return
                "c2_function.hh 490 2012-04-10 19:05:40Z marcus ";
        }

    /// \brief get versioning information for the source file
    /// \return the CVS Id string
        const std::string cvs_file_vers() const ;

public:
    /// \brief destructor
        virtual ~c2_function() {
                if(sampling_grid && !no_overwrite_grid) delete sampling_grid;
                if(root_info) delete root_info;
                if(owner_count) {
                        std::ostringstream outstr;
                        outstr << "attempt to delete an object with non-zero ownership in class ";
                        outstr << typeid(*this).name() << std::endl;
                        //throw c2_exception(outstr.str().c_str());
                }
        }


        virtual float_type value_with_derivatives(float_type x, float_type *yprime, float_type *yprime2) const// throw(c2_exception) =0 ;
         { return 0; };

        inline float_type operator () (float_type x) const// throw(c2_exception)
        { return value_with_derivatives(x, (float_type *)0, (float_type *)0); }

    inline float_type operator () (float_type x, float_type *yprime, float_type *yprime2) const// throw(c2_exception)
        { return value_with_derivatives(x, yprime, yprime2); }
        float_type find_root(float_type lower_bracket, float_type upper_bracket, float_type start,
        float_type value, int *error=0,
        float_type *final_yprime=0, float_type *final_yprime2=0 )const;// throw(c2_exception) ; // solve f(x)=value
        float_type partial_integrals(std::vector<float_type> xgrid, std::vector<float_type> *partials = 0,
          float_type abs_tol=1e-12, float_type rel_tol=1e-12, int derivs=2, bool adapt=true, bool extrapolate=true)
                const;// throw(c2_exception);

    /// \return sum of partial integrals, which is the definite integral from the first value in \a xgrid to the last.
        float_type integral(float_type amin, float_type amax, std::vector<float_type> *partials = 0,
             float_type abs_tol=1e-12, float_type rel_tol=1e-12, int derivs=2, bool adapt=true, bool extrapolate=true)
                const;// throw(c2_exception);

        c2_piecewise_function_p<float_type> *adaptively_sample(float_type amin, float_type amax,
                 float_type abs_tol=1e-12, float_type rel_tol=1e-12,
                 int derivs=2, std::vector<float_type> *xvals=0, std::vector<float_type> *yvals=0) const;// throw(c2_exception);

        /// \brief return the lower bound of the domain for this function as set by set_domain()
        inline float_type xmin() const { return fXMin; }
        /// \brief return the upper bound of the domain for this function as set by set_domain()
        inline float_type xmax() const { return fXMax; }
        /// \brief set the domain for this function.
        void set_domain(float_type amin, float_type amax) { fXMin=amin; fXMax=amax; }

        /// \brief this is a counter owned by the function but which can be used to monitor efficiency of algorithms.
        ///
        /// It is not maintained automatically in general!  The root finder, integrator, and sampler do increment it.
        /// \return number of evaluations logged since last reset.
        size_t get_evaluations() const { return evaluations; }
        /// \brief reset the counter
        void reset_evaluations()  const { evaluations=0; } // evaluations are 'invisible' to constant
        /// \brief count evaluations
        inline void increment_evaluations() const { evaluations++; }

        bool check_monotonicity(const std::vector<float_type> &data, const char message[])const;// throw(c2_exception);

        /// \param grid a vector of abscissas.  The contents is copied into an internal vector, so the \a grid can be discarded after passingin.
        virtual void set_sampling_grid(const std::vector<float_type> &grid); //throw(c2_exception);

        /// \brief get the sampling grid, which may be a null pointer
        /// \return pointer to the sampling grid
        std::vector<float_type> *get_sampling_grid_pointer() const { return sampling_grid; }

        virtual void get_sampling_grid(float_type amin, float_type amax, std::vector<float_type> &grid) const ;

        void preen_sampling_grid(std::vector<float_type> *result) const;

        void refine_sampling_grid(std::vector<float_type> &grid, size_t refinement) const;

        c2_function<float_type> &normalized_function(float_type amin, float_type amax, float_type norm=1.0)const;// throw(c2_exception);

        c2_function<float_type> &square_normalized_function(float_type amin, float_type amax, float_type norm=1.0)const;// throw(c2_exception);

        c2_function<float_type> &square_normalized_function(
                float_type amin, float_type amax, const c2_function<float_type> &weight, float_type norm=1.0)
                const;// throw(c2_exception);

        c2_sum_p<float_type> &operator + (const c2_function<float_type> &rhs)  const
                { return *new c2_sum_p<float_type>(*this, rhs); }

        c2_diff_p<float_type> &operator - (const c2_function<float_type> &rhs)  const
                { return *new c2_diff_p<float_type>(*this, rhs); }
        /// \brief factory function to create a c2_product_p from a regular algebraic expression.
        /// \param rhs the right-hand term of the product
        /// \return a new c2_function
        c2_product_p<float_type> &operator * (const c2_function<float_type> &rhs) const
                { return *new c2_product_p<float_type>(*this, rhs); }
        /// \brief factory function to create a c2_ratio_p from a regular algebraic expression.
        /// \param rhs the right-hand term of the ratio (the denominator)
        /// \return a new c2_function
        c2_ratio_p<float_type> &operator / (const c2_function<float_type> &rhs)  const
                { return *new c2_ratio_p<float_type>(*this, rhs); }
    /// \brief compose this function outside another.
    /// \param inner the inner function
    /// \return the composed function
        /// \anchor compose_operator
        c2_composed_function_p<float_type> & operator ()(const c2_function<float_type> &inner) const
                { return *new c2_composed_function_p<float_type>((*this), inner); }

        float_type get_trouble_point() const { return bad_x_point; }

        /// \brief increment our reference count.  Destruction is only legal if the count is zero.
        void claim_ownership() const { owner_count++; }
        /// \brief decrement our reference count. Do not destroy at zero.
        /// \return final owner count, to check whether object should disappear.
        size_t release_ownership_for_return() const //throw(c2_exception)
        {
                if(!owner_count) {
                        std::ostringstream outstr;
                        outstr << "attempt to release ownership of an unowned function in class ";
                        outstr << typeid(*this).name() << std::endl;
                        throw c2_exception(outstr.str().c_str());
                }
                owner_count--;
                return owner_count;
        }
        /// \brief decrement our reference count. If the count reaches zero, destroy ourself.
        void release_ownership() const// throw(c2_exception)
        {
                if(!release_ownership_for_return()) delete this;
        }
        /// \brief get the reference count, mostly for debugging
        /// \return the count
        size_t count_owners() const { return owner_count; }

protected:
        c2_function(const c2_function<float_type>  &src) : sampling_grid(0),
                no_overwrite_grid(false),
        fXMin(src.fXMin), fXMax(src.fXMax), root_info(0), owner_count(0)
        {} // copy constructor only copies domain, and is only for internal use
        c2_function() :
                        sampling_grid(0), no_overwrite_grid(0),
        fXMin(-std::numeric_limits<float_type>::max()),
        fXMax(std::numeric_limits<float_type>::max()), root_info(0), owner_count(0)
                {} // prevent accidental naked construction (impossible any since this has pure virtual methods)

        // this should only be called very early on, by a constructor, before anyone else
        // sets a sampling grid, or it will leak memory
        virtual void set_sampling_grid_pointer(std::vector<float_type> &grid)
                {
                        if (sampling_grid && !no_overwrite_grid) delete sampling_grid; // grid was ours, lose it.
                        sampling_grid=&grid; no_overwrite_grid=1;
                }

        std::vector<float_type> * sampling_grid;
        bool no_overwrite_grid;

        float_type fXMin, fXMax;
        mutable size_t evaluations;
        /// \brief this point may be used to record where a calculation ran into trouble
        mutable float_type bad_x_point;
public:
        /// \brief fill in a c2_fblock<float_type>... a shortcut for the integrator & sampler
        /// \param [in,out] fb the block to fill in with information
        inline void fill_fblock(c2_fblock<float_type> &fb) const// throw(c2_exception)
        {
                fb.y=value_with_derivatives(fb.x, &fb.yp, &fb.ypp);
                fb.ypbad=c2_isnan(fb.yp) || !c2_isfinite(fb.yp);
                fb.yppbad=c2_isnan(fb.ypp) || !c2_isfinite(fb.ypp);
                increment_evaluations();
        }

private:
        /// \brief the data element for the internal recursion stack for the sampler and integrator
        struct recur_item {
                c2_fblock<float_type> f1; size_t depth;
                float_type previous_estimate, abs_tol, step_sum;
                bool done;
                size_t f0index, f2index;
        };

        struct c2_integrate_recur {
                c2_fblock<float_type> *f0, *f1;
                float_type abs_tol, rel_tol, eps_scale, extrap_coef, extrap2, dx_tolerance, abs_tol_min;
                std::vector< recur_item > *rb_stack;
                int  derivs;
                bool adapt, extrapolate, inited;
        };

        /// \brief structure used to pass information recursively in sampler.
        ///
        struct c2_sample_recur {
                c2_fblock<float_type> *f0, *f1;
                float_type abs_tol, rel_tol, dx_tolerance, abs_tol_min;
                int derivs;
                c2_piecewise_function_p<float_type> *out;
                std::vector<float_type> *xvals, *yvals;
                std::vector< recur_item > *rb_stack;
                bool inited;
        };

        /// \brief structure used to hold root bracketing information
        ///
        struct c2_root_info {
                c2_fblock<float_type> lower, upper;
                bool inited;
        };


        float_type integrate_step(struct c2_integrate_recur &rb)const;// throw(c2_exception);
        void sample_step(struct c2_sample_recur &rb) const;// throw(c2_exception);

    mutable struct  c2_root_info *root_info;

        mutable size_t owner_count;
};

template <typename float_type=double> class c2_classic_function_p : public c2_function<float_type> {
public:
        /// \brief construct the container
        /// \param c_func a pointer to a conventional c-style function
        c2_classic_function_p(const float_type (*c_func)(float_type)) : c2_function<float_type>(), func(c_func)  {}

        /// \copydoc c2_function::value_with_derivatives
        /// Uses the internal function pointer set by set_function().
        virtual float_type value_with_derivatives(float_type x, float_type *yprime, float_type *yprime2) const// throw(c2_exception)
        {
                if(!func) throw c2_exception("c2_classic_function called with null function");
                if(yprime) *yprime=0;
                if(yprime2) *yprime2=0;
                return func(x);
        }
        virtual ~c2_classic_function_p() { }

protected:
        /// \brief pointer to our function
        const float_type (*func)(float_type);
};

template <typename float_type> class c2_const_ptr {
public:
        /// \brief construct the container with no function
        c2_const_ptr() : func(0)  {}
        /// \brief construct the container with a pre-defined function
        /// \param f the function to store
        c2_const_ptr(const c2_function<float_type> &f) : func(0)
                { this->set_function(&f); }
        /// \brief copy constructor
        /// \param src the container to copy
        c2_const_ptr(const c2_const_ptr<float_type> &src) : func(0)
                { this->set_function(src.get_ptr()); }
        /// \brief fill the container with a new function, or clear it with a null pointer
        /// \param f the function to store, releasing any previously held function
        void set_function(const c2_function<float_type> *f)
        {
                if(func) func->release_ownership();
                func=f;
                if(func) func->claim_ownership();
        }

        /// \brief fill the container from another container
        /// \param f the container to copy
        const c2_const_ptr<float_type> & operator =(const c2_const_ptr<float_type> &f)
                { this->set_function(f.get_ptr()); return f; }
        /// \brief fill the container with a function
        /// \param f the function
        const c2_function<float_type> & operator =(const c2_function<float_type> &f)
                { this->set_function(&f); return f; }

        void release_for_return() //throw(c2_exception)
                {
                        if(func) func->release_ownership_for_return();
                        func=0;
                }
        /// \brief clear the function
        ///
        /// Any attempt to use this c2_plugin_function_p throws an exception if the saved function is cleared.
        void unset_function(void) { this->set_function(0);  }
        /// \brief destructor
        ~c2_const_ptr() { this->set_function(0); }

        /// \brief get a reference to our owned function
        inline const c2_function<float_type> &get() const //throw(c2_exception)
        {
                if(!func) throw c2_exception("c2_ptr accessed uninitialized");
                return *func;
        }
        /// \brief get an unchecked pointer to our owned function
        inline const c2_function<float_type> *get_ptr() const { return func; }
        /// \brief get a checked pointer to our owned function
        inline const c2_function<float_type> *operator -> () const
                { return &get(); }
        /// \brief check if we have a valid function
        bool valid() const { return func != 0; }

        /// \brief type coercion operator which lets us use a pointer as if it were a const c2_function
        operator const c2_function<float_type>& () const { return this->get(); }

        float_type operator()(float_type x) const// throw(c2_exception)
        { return get()(x); }

        float_type operator()(float_type x, float_type *yprime, float_type *yprime2) const// throw(c2_exception)
                { return get().value_with_derivatives(x, yprime, yprime2); }

        c2_sum_p<float_type> &operator + (const c2_function<float_type> &rhs)  const// throw(c2_exception)
                { return *new c2_sum_p<float_type>(get(), rhs); }

        c2_diff_p<float_type> &operator - (const c2_function<float_type> &rhs)  const// throw(c2_exception)
                { return *new c2_diff_p<float_type>(get(), rhs); }

        c2_product_p<float_type> &operator * (const c2_function<float_type> &rhs) const// throw(c2_exception)
                { return *new c2_product_p<float_type>(get(), rhs); }

        c2_ratio_p<float_type> &operator / (const c2_function<float_type> &rhs)  const// throw(c2_exception)
                { return *new c2_ratio_p<float_type>(get(), rhs); }

        c2_composed_function_p<float_type> & operator ()(const c2_function<float_type> &inner) const// throw(c2_exception)
                { return *new c2_composed_function_p<float_type>(get(), inner); }

protected:
        const c2_function<float_type> * func;
};

template <typename float_type> class c2_ptr : public c2_const_ptr<float_type >
{
public:
        /// \brief construct the container with no function
        c2_ptr() : c2_const_ptr<float_type>()  {}
        /// \brief construct the container with a pre-defined function
        /// \param f the function to store
        c2_ptr(c2_function<float_type> &f) :
                c2_const_ptr<float_type>() { this->set_function(&f); }
        /// \brief copy constructor
        /// \param src the container to copy
        c2_ptr(const c2_ptr<float_type> &src) :
                c2_const_ptr<float_type>() { this->set_function(src.get_ptr()); }
        /// \brief get a checked pointer to our owned function
        inline c2_function<float_type> &get() const// throw(c2_exception)
                { return *const_cast<c2_function<float_type>*>(&c2_const_ptr<float_type>::get()); }
        /// \brief get an unchecked pointer to our owned function
        inline c2_function<float_type> *get_ptr() const
                { return const_cast<c2_function<float_type>*>(this->func); }
        /// \brief get a checked pointer to our owned function
        inline c2_function<float_type> *operator -> () const
                { return &get(); }
        /// \brief fill the container from another container
        /// \param f the container to copy
        const c2_ptr<float_type> & operator =(const c2_ptr<float_type> &f)
                { this->set_function(f.get_ptr()); return f; }
        /// \brief fill the container with a function
        /// \param f the function
        c2_function<float_type> & operator =(c2_function<float_type> &f)
                { this->set_function(&f); return f; }
private:
        /// \brief hidden non-const-safe version of operator=
        void operator =(const c2_const_ptr<float_type> &) { }
        /// \brief hidden non-const-safe version of operator=
        void operator =(const c2_function<float_type> &) { }
};

template <typename float_type, template <typename> class c2_class > class c2_typed_ptr : public c2_const_ptr<float_type> {
public:
        /// \brief construct the container with no function
        c2_typed_ptr() : c2_ptr<float_type>()  {}
        /// \brief construct the container with a pre-defined function
        /// \param f the function to store
        c2_typed_ptr(c2_class<float_type> &f)
                : c2_const_ptr<float_type>() { this->set_function(&f); }
        /// \brief copy constructor
        /// \param src the container to copy
        c2_typed_ptr(const c2_typed_ptr<float_type, c2_class> &src)
                : c2_const_ptr<float_type>() { this->set_function(src.get_ptr()); }

        /// \brief get a reference to our owned function
        inline c2_class<float_type> &get() const// throw(c2_exception)
                {
                        return *static_cast<c2_class<float_type> *>(const_cast<c2_function<float_type>*>(&c2_const_ptr<float_type>::get()));
                }
        /// \brief get a checked pointer to our owned function
        inline c2_class<float_type> *operator -> () const
                { return &get(); }
        /// \brief get an unchecked pointer to our owned function
        inline c2_class<float_type> *get_ptr() const
                { return static_cast<c2_class<float_type> *>(const_cast<c2_function<float_type>*>(this->func)); }
        /// \brief type coercion operator which lets us use a pointer as if it were a c2_function
        operator c2_class<float_type>&() const { return get(); }
        /// \brief fill the container from another container
        /// \param f the container to copy
        void operator =(const c2_typed_ptr<float_type, c2_class> &f)
                { this->set_function(f.get_ptr()); }
        /// \brief fill the container with a function
        /// \param f the function
        void operator =(c2_class<float_type> &f)
                { this->set_function(&f); }
private:
        /// \brief hidden downcasting version of operator=
        void operator =(const c2_const_ptr<float_type> &) { }
        /// \brief hidden downcasting version of operator=. Use an explicit dynamic_cast<c2_class<float_type>&>(f) if you need to try this.
        void operator =(const c2_function<float_type> &) { }
};

template <typename float_type=double> class c2_plugin_function_p :
        public c2_function<float_type> {
public:
        /// \brief construct the container with no function
        c2_plugin_function_p() : c2_function<float_type>(), func()  {}
        /// \brief construct the container with a pre-defined function
        c2_plugin_function_p(c2_function<float_type> &f) :
                c2_function<float_type>(),func(f)  { }
        /// \brief fill the container with a new function, or clear it with a null pointer
        /// and copy our domain
        void set_function(c2_function<float_type> *f)
                {
                        func.set_function(f);
                        if(f) this->set_domain(f->xmin(), f->xmax());
                }
        /// \copydoc c2_function::value_with_derivatives
        /// Uses the internal function pointer set by set_function().
        virtual float_type value_with_derivatives(float_type x, float_type *yprime, float_type *yprime2) const// throw(c2_exception)
                {
                        if(!func.valid()) throw c2_exception("c2_plugin_function_p called uninitialized");
                        return func->value_with_derivatives(x, yprime, yprime2);
                }
        /// \brief destructor
        virtual ~c2_plugin_function_p() { }

        /// \brief clear our function
        void unset_function() { func.unset_function(); }

        virtual void get_sampling_grid(float_type amin, float_type amax, std::vector<float_type> &grid) const
        {
                if(!func.valid()) throw c2_exception("c2_plugin_function_p called uninitialized");
                if(this->sampling_grid) c2_function<float_type>::get_sampling_grid(amin, amax, grid);
                else  func->get_sampling_grid(amin, amax, grid);
        }
protected:
        c2_ptr<float_type> func;
};

template <typename float_type=double> class c2_const_plugin_function_p : public c2_plugin_function_p<float_type> {
public:
        /// \brief construct the container with no function
        c2_const_plugin_function_p() : c2_plugin_function_p<float_type>()  {}
        /// \brief construct the container with a pre-defined function
        c2_const_plugin_function_p(const c2_function<float_type> &f) :
                c2_plugin_function_p<float_type>()  { this->set_function(&f); }
        /// \brief fill the container with a new function, or clear it with a null pointer
        void set_function(const c2_function<float_type> *f)
                { c2_plugin_function_p<float_type>::set_function(const_cast<c2_function<float_type>*>(f)); }
        /// \brief destructor
        virtual ~c2_const_plugin_function_p() { }

        /// \brief get a const reference to our owned function, for direct access
        const c2_function<float_type> &get() const// throw(c2_exception)
                { return this->func.get(); }
};

template <typename float_type=double> class c2_binary_function : public c2_function<float_type> {
public:
        ///  \brief function to manage the binary operation, used by c2_binary_function::value_with_derivatives()
    ///

        virtual float_type value_with_derivatives(float_type x, float_type *yprime, float_type *yprime2)// const throw (c2_exception)
    {
                if(stub) throw c2_exception("attempt to evaluate a c2_binary_function stub");
                return this->combine(*Left.get_ptr(), *Right.get_ptr(), x, yprime, yprime2);
    }

        /// \brief destructor releases ownership of member functions
        ///
        virtual ~c2_binary_function() { }

protected:

        c2_binary_function(
                        float_type (*combiner)(const c2_function<float_type> &left, const c2_function<float_type> &right,
                                   float_type x, float_type *yprime, float_type *yprime2),
                        const c2_function<float_type> &left,  const c2_function<float_type> &right) :
                        c2_function<float_type>(), combine(combiner), Left(left), Right(right), stub(false)
        {
                this->set_domain(
                                   (left.xmin() > right.xmin()) ? left.xmin() : right.xmin(),
                                   (left.xmax() < right.xmax()) ? left.xmax() : right.xmax()
                                   );
        }

    c2_binary_function(
           float_type (*combiner)(const c2_function<float_type> &left, const c2_function<float_type> &right,
                                                          float_type x, float_type *yprime, float_type *yprime2)
                ) : c2_function<float_type>(), combine(combiner), Left(), Right(), stub(true) { }

public:
        float_type (* const combine)(const c2_function<float_type> &left, const c2_function<float_type> &right,
                                                                 float_type x, float_type *yprime, float_type *yprime2);

protected:
        const c2_const_ptr<float_type> Left,  Right;
        /// \brief if true, we don't own any functions, we are just a source of a combining function.
        bool stub;

};

template <typename float_type=double> class c2_scaled_function_p : public c2_function<float_type> {
public:

        c2_scaled_function_p(const c2_function<float_type> &outer, float_type scale) :
                c2_function<float_type>(), func(outer), yscale(scale) { }

        /// \brief set a new scale factor
        /// \param scale the new factor
        void reset(float_type scale) { yscale=scale; }

        virtual float_type value_with_derivatives(float_type x, float_type *yprime, float_type *yprime2) const// throw (c2_exception)
    {
                float_type y=this->func->value_with_derivatives(x, yprime, yprime2);
                if(yprime) (*yprime)*=yscale;
                if(yprime2) (*yprime2)*=yscale;
                return y*yscale;
    }

protected:
    c2_scaled_function_p<float_type>() : func() {} // hide default constructor, since its use is almost always an error.
        /// \brief the scaling factor for the function
        const c2_const_ptr<float_type> func;
        float_type yscale;
};

template <typename float_type=double> class c2_cached_function_p : public c2_function<float_type> {
public:

        c2_cached_function_p(const c2_function<float_type> &f) : c2_function<float_type>(),
                func(f), init(false)  {}

        virtual float_type value_with_derivatives(float_type x, float_type *yprime, float_type *yprime2) const// throw(c2_exception)
    {
                if(!init || x != x0) {
                        y=this->func->value_with_derivatives(x, &yp, &ypp);
                        x0=x;
                        init=true;
                }
                if(yprime) *yprime=yp;
                if(yprime2) *yprime2=ypp;
                return y;
        }

protected:
    c2_cached_function_p() : func() {} // hide default constructor, since its use is almost always an error.
        const c2_const_ptr<float_type> func;
        mutable bool init;
        mutable float_type x0, y, yp, ypp;

};

template <typename float_type=double> class  c2_composed_function_p : public c2_binary_function<float_type> {
public:

        c2_composed_function_p(const c2_function<float_type> &outer, const c2_function<float_type> &inner) :
                c2_binary_function<float_type>(combine, outer, inner) { this->set_domain(inner.xmin(), inner.xmax()); }
        /// \brief Create a stub just for the combiner to avoid statics.
        c2_composed_function_p() : c2_binary_function<float_type>(combine) {}

        /// \brief execute math necessary to do composition
        static float_type combine(const c2_function<float_type> &left, const c2_function<float_type> &right,
                                                  float_type x, float_type *yprime, float_type *yprime2) //throw(c2_exception)
    {
                float_type y0, y1;
                if(yprime || yprime2) {
                        float_type yp0, ypp0, yp1, ypp1;
                        y0=right.value_with_derivatives(x, &yp0, &ypp0);
                        y1=left.value_with_derivatives(y0, &yp1, &ypp1);
                        if(yprime) *yprime=yp1*yp0;
                        if(yprime2) *yprime2=ypp0*yp1+yp0*yp0*ypp1;
                } else {
                        y0=right(x);
                        y1=left(y0);
                }
                return y1;
        }
};

template <typename float_type=double> class c2_sum_p : public c2_binary_function<float_type> {
public:

        c2_sum_p(const c2_function<float_type> &left, const c2_function<float_type> &right) : c2_binary_function<float_type>(combine, left, right) {}
        /// \brief Create a stub just for the combiner to avoid statics.
        c2_sum_p() : c2_binary_function<float_type>(combine) {} ; // create a stub just for the combiner to avoid statics

        /// \brief execute math necessary to do addition
        static float_type combine(const c2_function<float_type> &left, const c2_function<float_type> &right,
                                                  float_type x, float_type *yprime, float_type *yprime2) //throw(c2_exception)
    {
                float_type y0, y1;
                if(yprime || yprime2) {
                        float_type yp0, ypp0, yp1, ypp1;
                        y0=left.value_with_derivatives(x, &yp0, &ypp0);
                        y1=right.value_with_derivatives(x, &yp1, &ypp1);
                        if(yprime) *yprime=yp0+yp1;
                        if(yprime2) *yprime2=ypp0+ypp1;
                } else {
                        y0=left(x);
                        y1=right(x);
                }
                return y0+y1;
        }
};

template <typename float_type=double> class c2_diff_p : public c2_binary_function<float_type> {
public:

        c2_diff_p(const c2_function<float_type> &left, const c2_function<float_type> &right) : c2_binary_function<float_type>(combine, left, right) {}
        /// \brief Create a stub just for the combiner to avoid statics.
        c2_diff_p() : c2_binary_function<float_type>(combine) {} ; // create a stub just for the combiner to avoid statics

        /// \brief execute math necessary to do subtraction
        static float_type combine(const c2_function<float_type> &left, const c2_function<float_type> &right,
                                                          float_type x, float_type *yprime, float_type *yprime2) //throw(c2_exception)
    {
                float_type y0, y1;
                if(yprime || yprime2) {
                        float_type yp0, ypp0, yp1, ypp1;
                        y0=left.value_with_derivatives(x, &yp0, &ypp0);
                        y1=right.value_with_derivatives(x, &yp1, &ypp1);
                        if(yprime) *yprime=yp0-yp1;
                        if(yprime2) *yprime2=ypp0-ypp1;
                } else {
                        y0=left(x);
                        y1=right(x);
                }
                return y0-y1;
        }
};

template <typename float_type=double> class c2_product_p : public c2_binary_function<float_type> {
public:

        c2_product_p(const c2_function<float_type> &left, const c2_function<float_type> &right) : c2_binary_function<float_type>(combine, left, right) {}
        /// \brief Create a stub just for the combiner to avoid statics.
        c2_product_p() : c2_binary_function<float_type>(combine) {} ; // create a stub just for the combiner to avoid statics

        /// \brief execute math necessary to do multiplication
        static float_type combine(const c2_function<float_type> &left, const c2_function<float_type> &right,
                                                           float_type x, float_type *yprime, float_type *yprime2)  //throw(c2_exception)
    {
                float_type y0, y1;
                if(yprime || yprime2) {
                        float_type yp0, ypp0, yp1, ypp1;
                        y0=left.value_with_derivatives(x, &yp0, &ypp0);
                        y1=right.value_with_derivatives(x, &yp1, &ypp1);
                        if(yprime) *yprime=y1*yp0+y0*yp1;
                        if(yprime2) *yprime2=ypp0*y1+2.0*yp0*yp1+ypp1*y0;
                } else {
                        y0=left(x);
                        y1=right(x);
                }
                return y0*y1;
        }
};

template <typename float_type=double> class c2_ratio_p : public c2_binary_function<float_type> {
public:

        c2_ratio_p(const c2_function<float_type> &left, const c2_function<float_type> &right) : c2_binary_function<float_type>(combine, left, right) {}
        /// \brief Create a stub just for the combiner to avoid statics.
        c2_ratio_p() : c2_binary_function<float_type>(combine) {} ; // create a stub just for the combiner to avoid statics

        /// \brief execute math necessary to do division
        static float_type combine(const c2_function<float_type> &left, const c2_function<float_type> &right,
                                                          float_type x, float_type *yprime, float_type *yprime2) //throw(c2_exception)
    {
                float_type y0, y1;
                if(yprime || yprime2) {
                        float_type yp0, ypp0, yp1, ypp1;
                        y0=left.value_with_derivatives(x, &yp0, &ypp0);
                        y1=right.value_with_derivatives(x, &yp1, &ypp1);
                        if(yprime) *yprime=(yp0*y1-y0*yp1)/(y1*y1); // first deriv of ratio
                        if(yprime2) *yprime2=(y1*y1*ypp0+y0*(2*yp1*yp1-y1*ypp1)-2*y1*yp0*yp1)/(y1*y1*y1);
                } else {
                        y0=left(x);
                        y1=right(x);
                }
                return y0/y1;
        }

};

template <typename float_type> class c2_constant_p : public c2_function<float_type> {
public:
        c2_constant_p(float_type x) : c2_function<float_type>(), value(x) {}
        void reset(float_type val) { value=val; }
        virtual float_type value_with_derivatives(float_type, float_type *yprime, float_type *yprime2) const// throw(c2_exception)
        { if(yprime) *yprime=0; if(yprime2) *yprime2=0; return value; }

private:
                float_type value;
};

/// \brief a transformation of a coordinate, including an inverse
/// \ingroup transforms
template <typename float_type> class c2_transformation {
public:

        c2_transformation(bool transformed,
                   float_type (*xin)(float_type), float_type (*xinp)(float_type), float_type (*xinpp)(float_type), float_type (*xout)(float_type)
           ) :
                        fTransformed(transformed), fHasStaticTransforms(true),
                        pIn(xin), pInPrime(xinp), pInDPrime(xinpp), pOut(xout) { }

        c2_transformation(bool transformed) :
                        fTransformed(transformed), fHasStaticTransforms(false),
                        pIn(report_error), pInPrime(report_error), pInDPrime(report_error), pOut(report_error) { }
        /// \brief the destructor
        virtual ~c2_transformation() { }
        /// \brief flag to indicate if this transform is not the identity
        const bool fTransformed;
        /// \brief flag to indicate if the static function pointers can be used for efficiency
        const bool fHasStaticTransforms;

        float_type (* const pIn)(float_type);
        /// \brief non-virtual pointer to input X transform derivative
        float_type (* const pInPrime)(float_type);
        /// \brief non-virtual pointer to input X transform second derivative
        float_type (* const pInDPrime)(float_type);
        /// \brief non-virtual pointer to output X transform
        float_type (* const pOut)(float_type);

        /// \brief virtual input X transform
        virtual float_type fIn(float_type x) const { return pIn(x); }
        /// \brief virtual input X transform derivative
        virtual float_type fInPrime(float_type x) const { return pInPrime(x); }
        /// \brief virtual input X transform second derivative
        virtual float_type fInDPrime(float_type x) const { return pInDPrime(x); }
        /// \brief virtual output X transform
        virtual float_type fOut(float_type x) const { return pOut(x); }

protected:
        /// \brief utility function for unimplemented conversion
        static float_type report_error(float_type x)  { throw c2_exception("use of improperly constructed axis transform"); return x; }
        /// \brief utility function f(x)=x useful in axis transforms
        static float_type ident(float_type x)  { return x; }
        /// \brief utility function f(x)=1 useful in axis transforms
        static float_type one(float_type)  { return 1; }
        /// \brief utility function f(x)=0 useful in axis transforms
        static float_type zero(float_type)  { return 0; }
        /// \brief utility function f(x)=1/x useful in axis transforms
        static float_type recip(float_type x)  { return 1.0/x; }
        /// \brief utility function f(x)=-1/x**2 useful in axis transforms
        static float_type recip_prime(float_type x)  { return -1/(x*x); }
        /// \brief utility function f(x)=2/x**3 useful in axis transforms
        static float_type recip_prime2(float_type x)  { return 2/(x*x*x); }

};

/// \brief the identity transform
/// \ingroup transforms
template <typename float_type> class c2_transformation_linear : public c2_transformation<float_type> {
public:
        /// \brief constructor
        c2_transformation_linear() : c2_transformation<float_type>(false, this->ident, this->one, this->zero, this->ident) { }
        /// \brief destructor
        virtual ~c2_transformation_linear() { }
};
/// \brief log axis transform
/// \ingroup transforms
template <typename float_type> class c2_transformation_log : public c2_transformation<float_type> {
public:
        /// \brief constructor
        c2_transformation_log() : c2_transformation<float_type>(true, std::log, this->recip, this->recip_prime, std::exp) { }
        /// \brief destructor
        virtual ~c2_transformation_log() { }
};
/// \brief reciprocal axis transform
/// \ingroup transforms
template <typename float_type> class c2_transformation_recip : public c2_transformation<float_type> {
public:
        /// \brief constructor
        c2_transformation_recip() : c2_transformation<float_type>(true, this->recip, this->recip_prime, this->recip_prime2, this->recip) { }
        /// \brief destructor
        virtual ~c2_transformation_recip() { }
};

template <typename float_type>
        class c2_function_transformation {
public:
        /// \brief construct this from two c2_transformation instances
        /// \param xx the X axis transform
        /// \param yy the Y axis transform
        c2_function_transformation(
                const c2_transformation<float_type> &xx, const c2_transformation<float_type> &yy) :
                isIdentity(!(xx.fTransformed || yy.fTransformed)), X(xx), Y(yy) { }
        /// \brief destructor
        virtual ~c2_function_transformation() { delete &X; delete &Y; }

        virtual float_type evaluate(float_type xraw,
                float_type y, float_type yp0, float_type ypp0,
                float_type *yprime, float_type *yprime2) const;
        /// \brief flag indicating of the transform is the identity, and can be skipped for efficiency
        const bool isIdentity;
        /// \brief the X axis transform
        const c2_transformation<float_type> &X;
        /// \brief the Y axis transform
        const c2_transformation<float_type> &Y;
};

/// \brief a transformation of a function in and out of lin-lin space
///
/// \ingroup transforms
template <typename float_type> class c2_lin_lin_function_transformation :
        public c2_function_transformation<float_type> {
public:
        c2_lin_lin_function_transformation() :
                c2_function_transformation<float_type>(
                        *new c2_transformation_linear<float_type>,
                        *new c2_transformation_linear<float_type>
                ) { }
        virtual ~c2_lin_lin_function_transformation() { }
};

/// \brief a transformation of a function in and out of log-log space
///
/// \ingroup transforms
template <typename float_type> class c2_log_log_function_transformation :
        public c2_function_transformation<float_type> {
public:
        c2_log_log_function_transformation() :
                c2_function_transformation<float_type>(
                        *new c2_transformation_log<float_type>,
                        *new c2_transformation_log<float_type>
                ) { }
        virtual ~c2_log_log_function_transformation() { }
};

/// \brief a transformation of a function in and out of lin-log space
///
/// \ingroup transforms
template <typename float_type> class c2_lin_log_function_transformation :
        public c2_function_transformation<float_type> {
public:
        c2_lin_log_function_transformation() :
                c2_function_transformation<float_type>(
                        *new c2_transformation_linear<float_type>,
                        *new c2_transformation_log<float_type>
                ) { }
        virtual ~c2_lin_log_function_transformation() { }
};

/// \brief a transformation of a function in and out of log-lin space
///
/// \ingroup transforms
template <typename float_type> class c2_log_lin_function_transformation :
        public c2_function_transformation<float_type> {
public:
        c2_log_lin_function_transformation() :
                c2_function_transformation<float_type>(
                        *new c2_transformation_log<float_type>,
                        *new c2_transformation_linear<float_type>
                ) { }
        virtual ~c2_log_lin_function_transformation() { }
};

/// \brief a transformation of a function in and out of Arrhenius (1/x vs. log(y)) space
///
/// \ingroup transforms
template <typename float_type> class c2_arrhenius_function_transformation :
        public c2_function_transformation<float_type> {
public:
        c2_arrhenius_function_transformation() :
                c2_function_transformation<float_type>(
                        *new c2_transformation_recip<float_type>,
                        *new c2_transformation_log<float_type>
                ) { }
        virtual ~c2_arrhenius_function_transformation() { }
};

/**
    \brief create a cubic spline interpolation of a set of (x,y) pairs
        \ingroup interpolators
    This is one of the main reasons for c2_function objects to exist.

    It provides support for cubic spline interpolation of data provides from tables of \a x, \a y pairs.
    It supports automatic, transparent linearization of the data before storing in its tables (through
    subclasses such as
    log_lin_interpolating_function, lin_log_interpolating_function, and
    log_log_interpolating_function) to permit very high accuracy representations of data which have a suitable
    structure.  It provides utility functions LinearInterpolatingGrid() and LogLogInterpolatingGrid()
    to create grids for mapping other functions onto a arithmetic or geometric grid.

    In its simplest form, an untransformed cubic spline of a data set, using natural boundary conditions
    (vanishing second derivative), is created as: \n
    \code
        c2_ptr<double> c2p;
        c2_factory<double> c2;
    std::vector<double> xvals(10), yvals(10);
    // < fill in xvals and yvals >
    c2p myfunc=c2.interpolating_function().load(xvals, yvals,true,0,true,0);
    // and it can be evaluated at a point for its value only by:
    double y=myfunc(x);
    // or it can be evaluated with its derivatives by
    double yprime, yprime2;
    double y=myfunc(x,&yprime, &yprime2);
    \endcode

        The factory function c2_factory::interpolating_function() creates *new interpolating_function_p()
*/

template <typename float_type=double> class interpolating_function_p  : public c2_function<float_type> {
public:
    /// \brief an empty linear-linear cubic-spline interpolating_function_p
    ///
    /// lots to say here, but see Numerical Recipes for a discussion of cubic splines.
        ///
        interpolating_function_p() : c2_function<float_type>(),
                fTransform(*new  c2_lin_lin_function_transformation<float_type>) { }

    /// \brief an empty cubic-spline interpolating_function_p with a specific transform
    ///
        interpolating_function_p(const c2_function_transformation<float_type> &transform) : c2_function<float_type>(),
                fTransform(transform) { }


        interpolating_function_p<float_type> & load(const std::vector<float_type> &x, const std::vector<float_type> &f,
                          bool lowerSlopeNatural, float_type lowerSlope,
                          bool upperSlopeNatural, float_type upperSlope, bool splined=true
        ) ;//throw(c2_exception);


        interpolating_function_p<float_type> & load_pairs(
                          std::vector<std::pair<float_type, float_type> > &data,
                          bool lowerSlopeNatural, float_type lowerSlope,
                          bool upperSlopeNatural, float_type upperSlope, bool splined=true
        ) ;//throw(c2_exception);


        interpolating_function_p<float_type> &  sample_function(const c2_function<float_type> &func,
                          float_type amin, float_type amax, float_type abs_tol, float_type rel_tol,
                          bool lowerSlopeNatural, float_type lowerSlope,
                          bool upperSlopeNatural, float_type upperSlope
                        ) ;//throw(c2_exception);


        interpolating_function_p<float_type> & load_random_generator_function(
                const std::vector<float_type> &bincenters, const c2_function<float_type> &binheights)
                ;//throw(c2_exception);

        interpolating_function_p<float_type> & load_random_generator_bins(
                const std::vector<float_type> &bins, const std::vector<float_type> &binheights, bool splined=true)
                ;//throw(c2_exception);

        virtual float_type value_with_derivatives(float_type x, float_type *yprime, float_type *yprime2) const;// throw(c2_exception);

        /// \brief destructor
        virtual ~interpolating_function_p() { delete &fTransform; }

        /// \brief create a new, empty interpolating function of this type (virtual constructor)
        virtual interpolating_function_p<float_type> &clone() const// throw(c2_exception)
        {        return *new interpolating_function_p<float_type>(); }

        void get_data(std::vector<float_type> &xvals, std::vector<float_type> &yvals) const;// throw() ;

        void get_internal_data(
        std::vector<float_type> &xvals,
        std::vector<float_type> &yvals,
        std::vector<float_type> &y2vals) const
        { xvals=X; yvals=F; y2vals=y2; }

        void set_lower_extrapolation(float_type bound);

        void set_upper_extrapolation(float_type bound);

        interpolating_function_p <float_type> & unary_operator(const c2_function<float_type> &source) const;

        interpolating_function_p <float_type> & binary_operator(const c2_function<float_type> &rhs,
           const c2_binary_function<float_type> *combining_stub
           ) const;

    interpolating_function_p <float_type> & add_pointwise (const c2_function<float_type> &rhs) const {
                return binary_operator(rhs, new c2_sum_p<float_type>()); }

        interpolating_function_p <float_type> & subtract_pointwise (const c2_function<float_type> &rhs) const {
                return binary_operator(rhs, new c2_diff_p<float_type>()); }

        interpolating_function_p <float_type> & multiply_pointwise (const c2_function<float_type> &rhs) const {
                return binary_operator(rhs, new c2_product_p<float_type>()); }

        interpolating_function_p <float_type> & divide_pointwise (const c2_function<float_type> &rhs) const {
                return binary_operator(rhs, new c2_ratio_p<float_type>()); }

        void clone_data(const interpolating_function_p <float_type> &rhs) {
                Xraw=rhs.Xraw; X=rhs.X; F=rhs.F; y2=rhs.y2;
                set_sampling_grid_pointer(Xraw);
        }

        const c2_function_transformation<float_type> &fTransform;

protected:
    /// \brief create the spline coefficients
        void spline(
                          bool lowerSlopeNatural, float_type lowerSlope,
                          bool upperSlopeNatural, float_type upperSlope
                        ) ;//throw(c2_exception);

        // This is for sorting the data. It must be static if it's going to be a class member.
        static bool comp_pair(std::pair<float_type,float_type> const &i, std::pair<float_type,float_type> const &j) {return i.first<j.first;}

    std::vector<float_type> Xraw, X, F, y2;
        c2_const_ptr<float_type> sampler_function;
        bool xInverted;
    mutable size_t lastKLow;
};

///        The factory function c2_factory::log_lin_interpolating_function() creates *new log_lin_interpolating_function_p()
template <typename float_type=double> class log_lin_interpolating_function_p : public interpolating_function_p <float_type> {
public:
    /// \brief an empty log-linear cubic-spline interpolating_function_p
    ///
        log_lin_interpolating_function_p() : interpolating_function_p<float_type>(*new c2_log_lin_function_transformation<float_type>)
                { }
        /// \brief create a new, empty interpolating function of this type (virtual constructor)
        virtual interpolating_function_p<float_type> &clone() const// throw(c2_exception)
                {        return *new log_lin_interpolating_function_p<float_type>(); }
};

///        The factory function c2_factory::lin_log_interpolating_function() creates *new lin_log_interpolating_function_p()
template <typename float_type=double> class lin_log_interpolating_function_p : public interpolating_function_p <float_type> {
public:
    /// \brief an empty linear-log cubic-spline interpolating_function_p
    ///
        lin_log_interpolating_function_p() : interpolating_function_p<float_type>(*new c2_lin_log_function_transformation<float_type>)
                { }
        /// \brief create a new, empty interpolating function of this type (virtual constructor)
        virtual interpolating_function_p<float_type> &clone() const// throw(c2_exception)
                {        return *new lin_log_interpolating_function_p<float_type>(); }
};



///        The factory function c2_factory::log_log_interpolating_function() creates *new log_log_interpolating_function_p()
template <typename float_type=double> class log_log_interpolating_function_p : public interpolating_function_p <float_type> {
public:
    /// \brief an empty log-log cubic-spline interpolating_function_p
    ///
        log_log_interpolating_function_p() : interpolating_function_p<float_type>(*new c2_log_log_function_transformation<float_type>)
                { }
        /// \brief create a new, empty interpolating function of this type (virtual constructor)
        virtual interpolating_function_p<float_type> &clone() const// throw(c2_exception)
        {        return *new log_log_interpolating_function_p<float_type>(); }
};

template <typename float_type=double> class arrhenius_interpolating_function_p : public interpolating_function_p <float_type> {
public:
    /// \brief an empty arrhenius cubic-spline interpolating_function_p
    ///
        arrhenius_interpolating_function_p() : interpolating_function_p<float_type>(*new c2_arrhenius_function_transformation<float_type>)
                { }
        /// \brief create a new, empty interpolating function of this type (virtual constructor)
        virtual interpolating_function_p<float_type> &clone() const// throw(c2_exception)
        {        return *new arrhenius_interpolating_function_p<float_type>(); }
};


/// The factory function c2_factory::sin() creates *new c2_sin_p
template <typename float_type=double> class c2_sin_p : public c2_function<float_type> {
public:
        /// \brief constructor.
        c2_sin_p() : c2_function<float_type>() {}

        virtual float_type value_with_derivatives(float_type x, float_type *yprime, float_type *yprime2) const// throw(c2_exception)
        { float_type q=std::sin(x); if(yprime) *yprime=std::cos(x); if(yprime2) *yprime2=-q; return q; }

        virtual void get_sampling_grid(float_type amin, float_type amax,  std::vector<float_type> &grid) const;
};

template <typename float_type=double> class c2_cos_p : public c2_sin_p<float_type> {
public:
        /// \brief constructor.
        c2_cos_p() : c2_sin_p<float_type>() {}

        virtual float_type value_with_derivatives(float_type x, float_type *yprime, float_type *yprime2) const// throw(c2_exception)
        { float_type q=std::cos(x); if(yprime) *yprime=-std::sin(x); if(yprime2) *yprime2=-q; return q; }
};

/// \brief compute tan(x) with its derivatives.
/// \ingroup math_functions
///
/// The factory function c2_factory::tan() creates *new c2_tan_p
template <typename float_type=double> class c2_tan_p : public c2_function<float_type> {
public:
        /// \brief constructor.
        c2_tan_p() : c2_function<float_type>() {}

        virtual float_type value_with_derivatives(float_type x, float_type *yprime, float_type *yprime2) const// throw(c2_exception)
        {
                float_type c=std::cos(x), ss=std::sin(x);
                float_type t=ss/c;
                float_type yp=1/(c*c);
                if(yprime) {*yprime=yp;}
                if(yprime2) {*yprime2=2*t*yp;}
                return t;
        }
};

/// \brief compute log(x) with its derivatives.
/// \ingroup math_functions
///
/// The factory function c2_factory::log() creates *new c2_log_p
template <typename float_type=double> class c2_log_p : public c2_function<float_type> {
public:
        /// \brief constructor.
        c2_log_p() : c2_function<float_type>() {}

        virtual float_type value_with_derivatives(float_type x, float_type *yprime, float_type *yprime2) const// throw(c2_exception)
        { if(yprime) *yprime=1.0/x; if(yprime2) *yprime2=-1.0/(x*x); return std::log(x); }
};

/// \brief compute exp(x) with its derivatives.
/// \ingroup math_functions
///
/// The factory function c2_factory::exp() creates *new c2_exp_p
template <typename float_type=double>  class c2_exp_p : public c2_function<float_type> {
public:
        /// \brief constructor.
        c2_exp_p() : c2_function<float_type>() {}

        virtual float_type value_with_derivatives(float_type x, float_type *yprime, float_type *yprime2) const// throw(c2_exception)
        { float_type q=std::exp(x); if(yprime) *yprime=q; if(yprime2) *yprime2=q; return q; }
};

/// \brief compute sqrt(x) with its derivatives.
/// \ingroup math_functions
///
/// The factory function c2_factory::sqrt() creates *new c2_sqrt_p()
template <typename float_type=double> class c2_sqrt_p : public c2_function<float_type> {
public:
        /// \brief constructor.
        c2_sqrt_p() : c2_function<float_type>() {}

        virtual float_type value_with_derivatives(float_type x, float_type *yprime, float_type *yprime2) const// throw(c2_exception)
        { float_type q=std::sqrt(x); if(yprime) *yprime=0.5/q; if(yprime2) *yprime2=-0.25/(x*q); return q; }
};

/// \brief compute scale/x with its derivatives.
/// \ingroup parametric_functions
///
/// The factory function c2_factory::recip() creates *new c2_recip_p
template <typename float_type=double> class c2_recip_p : public c2_function<float_type> {
public:
        /// \brief constructor.
        c2_recip_p(float_type scale) : c2_function<float_type>(), rscale(scale) {}

        virtual float_type value_with_derivatives(float_type x, float_type *yprime, float_type *yprime2) const// throw(c2_exception)
        {
                float_type q=1.0/x;
                float_type y=rscale*q;
                if(yprime) *yprime=-y*q;
                if(yprime2) *yprime2=2*y*q*q;
                return y;
        }
        /// \brief reset the scale factor
        /// \param scale the new numerator
        void reset(float_type scale) { rscale=scale; }
private:
        float_type rscale;
};

/// \brief compute x with its derivatives.
/// \ingroup math_functions
///
/// The factory function c2_factory::identity() creates *new c2_identity_p
template <typename float_type=double> class c2_identity_p : public c2_function<float_type> {
public:
        /// \brief constructor.
        c2_identity_p() : c2_function<float_type>() {}

        virtual float_type value_with_derivatives(float_type x, float_type *yprime, float_type *yprime2) const// throw(c2_exception)
        { if(yprime) *yprime=1.0; if(yprime2) *yprime2=0; return x; }
};

template <typename float_type=double> class c2_linear_p : public c2_function<float_type> {
public:

        c2_linear_p(float_type x0, float_type y0, float_type slope) :
                c2_function<float_type>(), xint(x0), intercept(y0), m(slope) {}

        void reset(float_type x0, float_type y0, float_type slope) { xint=x0; intercept=y0; m=slope; }
        virtual float_type value_with_derivatives(float_type x, float_type *yprime, float_type *yprime2) const// throw(c2_exception)
        { if(yprime) *yprime=m; if(yprime2) *yprime2=0; return m*(x-xint)+intercept; }

private:
                float_type xint, intercept, m;
protected:
                c2_linear_p() {} // do not allow naked construction... it is usually an accident.
};


template <typename float_type=double> class c2_quadratic_p : public c2_function<float_type> {
public:

        c2_quadratic_p(float_type x0, float_type y0, float_type xcoef, float_type x2coef) :
                c2_function<float_type>(), intercept(y0), center(x0), a(x2coef), b(xcoef) {}

        void reset(float_type x0, float_type y0, float_type xcoef, float_type x2coef) { intercept=y0; center=x0; a=x2coef; b=xcoef; }
        virtual float_type value_with_derivatives(float_type x, float_type *yprime, float_type *yprime2) const// throw(c2_exception)
        { float_type dx=x-center; if(yprime) *yprime=2*a*dx+b; if(yprime2) *yprime2=2*a; return a*dx*dx+b*dx+intercept; }

private:
                float_type intercept, center, a, b;
protected:
                c2_quadratic_p() {} // do not allow naked construction... it is usually an accident.
};

template <typename float_type=double> class c2_power_law_p : public c2_function<float_type> {
public:
    /// \brief Construct the operator
    /// \param scale the multipler
    /// \param power the exponent
        c2_power_law_p(float_type scale, float_type power) :
                c2_function<float_type>(), a(scale), b(power) {}
    /// \brief Modify the mapping after construction
    /// \param scale the new multipler
    /// \param power the new exponent
        void reset(float_type scale, float_type power) { a=scale; b=power; }
        virtual float_type value_with_derivatives(float_type x, float_type *yprime, float_type *yprime2) const// throw(c2_exception)
        { float_type q=a*std::pow(x,b-2); if(yprime) *yprime=b*q*x; if(yprime2) *yprime2=b*(b-1)*q; return q*x*x; }

private:
                float_type a, b;
protected:
                c2_power_law_p() {} // do not allow naked construction... it is usually an accident.
};


template <typename float_type=double> class c2_inverse_function_p : public c2_function<float_type> {
public:
    /// \brief Construct the operator
    /// \param source the function to be inverted
        c2_inverse_function_p(const c2_function<float_type> &source);
    virtual float_type value_with_derivatives(float_type x, float_type *yprime, float_type *yprime2) const;// throw(c2_exception);

    /// \brief give the function a hint as to where to look for its inverse
    /// \param hint the likely value of the inverse, which defaults to whatever the evaluation returned.
    void set_start_hint(float_type hint) const { start_hint=hint; }

    virtual float_type get_start_hint(float_type x) const
                { return hinting_function.valid()? hinting_function(x) : start_hint; }

        void set_hinting_function(const c2_function<float_type> *hint_func)
                { hinting_function.set_function(hint_func); }

        void set_hinting_function(const c2_const_ptr<float_type> hint_func)
                { hinting_function=hint_func; }

protected:
        c2_inverse_function_p() {} // do not allow naked construction... it is usually an accident.
        mutable float_type start_hint;
        const c2_const_ptr<float_type> func;
        c2_const_ptr<float_type> hinting_function;
};


template <typename float_type=double>  class accumulated_histogram : public interpolating_function_p <float_type> {
public:

        accumulated_histogram(const std::vector<float_type>binedges, const std::vector<float_type> binheights,
                                                 bool normalize=false, bool inverse_function=false, bool drop_zeros=true);

};

template <typename float_type=double> class c2_connector_function_p : public c2_function<float_type> {
public:

        c2_connector_function_p(float_type x0, const c2_function<float_type> &f0, float_type x2, const c2_function<float_type> &f2,
                        bool auto_center, float_type y1);

        c2_connector_function_p(
                float_type x0, float_type y0, float_type yp0, float_type ypp0,
                float_type x2, float_type y2, float_type yp2, float_type ypp2,
                bool auto_center, float_type y1);

        c2_connector_function_p(
                const c2_fblock<float_type> &fb0,
                const c2_fblock<float_type> &fb2,
                bool auto_center, float_type y1);

        /// \brief destructor
        virtual ~c2_connector_function_p();
        virtual float_type value_with_derivatives(float_type x, float_type *yprime, float_type *yprime2) const;// throw (c2_exception);
protected:
        /// \brief fill container numerically
        void init(
                const c2_fblock<float_type> &fb0,
                const c2_fblock<float_type> &fb2,
                bool auto_center, float_type y1);

        float_type fhinv, fy1, fa, fb, fc, fd, fe, ff;
};

template <typename float_type=double> class c2_piecewise_function_p : public c2_function<float_type> {
public:
        /// \brief construct the container
        c2_piecewise_function_p();
        /// \brief destructor
        virtual ~c2_piecewise_function_p();
        virtual float_type value_with_derivatives(float_type x, float_type *yprime, float_type *yprime2) const;// throw (c2_exception);

        void append_function(const c2_function<float_type> &func) ;//throw (c2_exception);
protected:
        std::vector<c2_const_ptr<float_type> > functions;
        mutable int lastKLow;
};

#include "c2_function.icc"

#endif
