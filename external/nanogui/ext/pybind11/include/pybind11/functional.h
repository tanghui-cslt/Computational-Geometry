/*
    pybind11/functional.h: std::function<> support

    Copyright (c) 2016 Wenzel Jakob <wenzel.jakob@epfl.ch>

    All rights reserved. Use of this source code is governed by a
    BSD-style license that can be found in the LICENSE file.
*/

#pragma once

#include "pybind11.h"
#include <functional>

NAMESPACE_BEGIN(pybind11)
NAMESPACE_BEGIN(detail)

template <typename Return, typename... Args>
struct type_caster<std::function<Return(Args...)>> {
    using type = std::function<Return(Args...)>;
    using retval_type = conditional_t<std::is_same<Return, void>::value, void_type, Return>;
    using function_type = Return (*) (Args...);

public:
    bool load(handle src_, bool) {
        if (src_.is_none())
            return true;

        src_ = detail::get_function(src_);
        if (!src_ || !PyCallable_Check(src_.ptr()))
            return false;

        /*
           When passing a C++ function as an argument to another C++
           function via Python, every function call would normally involve
           a full C++ -> Python -> C++ roundtrip, which can be prohibitive.
           Here, we try to at least detect the case where the function is
           stateless (i.e. function pointer or lambda function without
           captured variables), in which case the roundtrip can be avoided.
         */
        if (PyCFunction_Check(src_.ptr())) {
            auto c = reinterpret_borrow<capsule>(PyCFunction_GET_SELF(src_.ptr()));
            auto rec = (function_record *) c;

            if (rec && rec->is_stateless && rec->data[1] == &typeid(function_type)) {
                struct capture { function_type f; };
                value = ((capture *) &rec->data)->f;
                return true;
            }
        }

        auto src = reinterpret_borrow<object>(src_);
        value = [src](Args... args) -> Return {
            gil_scoped_acquire acq;
            object retval(src(std::forward<Args>(args)...));
            /* Visual studio 2015 parser issue: need parentheses around this expression */
            return (retval.template cast<Return>());
        };
        return true;
    }

    template <typename Func>
    static handle cast(Func &&f_, return_value_policy policy, handle /* parent */) {
        if (!f_)
            return none().inc_ref();

        auto result = f_.template target<function_type>();
        if (result)
            return cpp_function(*result, policy).release();
        else
            return cpp_function(std::forward<Func>(f_), policy).release();
    }

    PYBIND11_TYPE_CASTER(type, _("Callable[[") +
            argument_loader<Args...>::arg_names() + _("], ") +
            make_caster<retval_type>::name() +
            _("]"));
};

NAMESPACE_END(detail)
NAMESPACE_END(pybind11)
