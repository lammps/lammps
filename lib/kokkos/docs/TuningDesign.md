#Kokkos Tuning

This is a design document describing the motivation, ideas, design, and prototype implementation of the Kokkos Tuning System

## Motivation

Currently, Kokkos makes a lot of decisions about tuning parameters (CUDA block sizes, different kernel implementations)
by picking an option that results in the best performance for the widest array of applications and architectures at the
time the choice is made. This approach leaves performance on the table, and appears increasingly untenable as the number
of architectures and applications grows, and as software versions change.

The Kokkos team would like to instead open up the ability to set the parameters as part of the tooling system so that
these parameters can be tuned for individual applications across all the architectures they might run on. In order to match the
feel of past Kokkos tooling efforts, we'd like to achieve this with a callback system.

## Ideas

A Kokkos Tuning system should be as small as is wise while achieving the following goals

1. Expose to tools enough data about the _context_ of the running application to tune intelligently. In autotuning terms, decribe the _features_
2. Expose to tools enough data about tuning parameters that they might know how to optimize what they're asked to
3. Expose to applications an interface that they might inform a tool about their current application context
4. Expose to tools the results of their choices
5. No perturbation of Kokkos Core when this system is disabled

Shared among the first three of these goals is a need for some way to describe the semantics of variables (tuning parameters, context variables)
internal to Kokkos or an application to an outside tool.

### Semantics of Variables

I think it's best to talk about the semantics of variables with concrete examples.

Suppose Kokkos wants a tool to choose a block size for it. Suppose all the application context is perfectly understood, that the tool knows
that the application has 10,000,000 particles active and that it's running a kernel called "make_particles_go," which is a parallel_for in
the "cuda" execution space. Even with this knowledge, the tool needs to know several things about what a block size _is_ for this to be generic and practical

1. Is it an integer value? A float? A string? (Type)
2. Relatedly, what are the mathematical semantics which are valid for it? Is it something
for which a list can be sorted? Do the distances between items in a sorted list make sense?
If I divide two values, does the ratio have some meaning? (semantics)
3. What are the valid choices for this value? Is a block size of -128 okay? How about 7? (candidates)

Semantics (as always) are likely the source of the most confusion here, so a bit of detail is good. Here I'm leaning heavily on the field
of statistics to enable tools to do intelligent searching. If ordering doesn't make sense, if a value is "categorical", the only thing
a tool can do is try all possible values for a tuning value. If they're ordered (ordinal), the search can take advantage of this by
using the concept of a directional search. If the distances between elements matter (interval data) you can cheat with things like
bisection. Finally if ratios matter you can play games where you increase by a factor of 10 in your searches. Note that one good point in favor of this design is that it matches up nicely with scikit-opt (a happy accident).

In describing the candidate values in (3), users have two options: sets or ranges. A set has a number of entries of the given type, a range has lower and upper bounds and a step size.

Claim: the combination of context, candidates, semantics, and types gives a tool enough to intelligently explore the search space of
tuning parameters

### Context

Suppose a tool perfectly understands what a block size is. To effectively tune one, it needs to know something about the application.

In a trivial case, the tool knows absolutely nothing other than candidate values for the block size, and tries to make a choice that optimizes across all
invocations of kernels. This isn't _that_ far from what Kokkos does now, so it's not unreasonable for this to produce decent results.
That said, we could quickly add some context from Kokkos, stuff like the name and type of the kernel, the execution space, all with the semantic information described above. That way a tuning tool could differentiate based on all the information available to Kokkos. Going a little further, we could expose this ability to provide context to our applications. What if the tools wasn't just tuning to the fact that the kernel name was "GEMM", but that "matrix_size" was a million? Or that "live_particles" had a
certain value? The more (relevant) context we provide to a tool, the better it will be able to tune.


### Intended Tool Workflow

Okay, so a tool knows what it's tuning, and it knows the context of the application well enough to do clever ML things, all of this with happy semantic information so that everything make . What should a workflow look like? A tool should

1) Listen to declarations about the semantics of context and tuning variables
2) Make tuning decisions
3) Measure their feedback
4) Get better at (2)

The easier we make this loop, the better

## Design

The design of this system is intended to reflect the above ideas with the minimal necessary additions to make the mechanics work. This section is almost entirely describing the small holes in the above descriptions. Variable declaration works exactly as described above, except that we associate types and associated values using a type_id with each type at declaration time.

Any time a value of a variable is declared (context) or requested (tuning), it is also associated with a context ID that says how long that declaration is valid for. So if a user sees

```c++
startContext(contextId(0))
declare_value("is_safe_to_push_button",true,contextId(0));
foo();
endContext(contextId(0));
bar();
```

They should know in `bar` that it is no longer safe to push the button. Similarly, if tools have provided tuning values to contextId(0), when contextId(0) ends, that is when the tool takes measurements related to those tuning values and learns things. *For most tools, when they see a call to startContext associated with a contextId, they'll do a starting measurement, and at endContext they'll stop that measurement*.

One ugly bit of semantic complexity is in variables with complicated sets of candidates. Taking the exmaple of GPU block size, for different kernels an application might have different sets of valid block sizes. This means that while "block size" might make sense as a type, there could be different types, "block_sizes_up_to_1024," "block_sizes_up_to_2048," that cover the concept of block size. In our experience every solution to this problem is ugly, our alternate answers were much uglier.

## Implementation

This section describes the implementation.

If you're writing a tool, you care about tool implementation.

If you want tools to know about information from your application, you care about application implementation

If you're a Kokkos developer, you care about the application implementation and Kokkos implementation

### Tool implementation

In the past, tools have responded to the [profiling hooks in Kokkos](https://github.com/kokkos/kokkos-tools/wiki/Profiling-Hooks). This effort adds to that, there are now a few more functions (note that I'm using the C names for types. In general you can replace Kokkos_Tools_ with Kokkos::Tools:: in C++ tools)


```c++
void kokkosp_declare_output_type(const char* name, const size_t id, Kokkos_Tools_VariableInfo& info);
```

Declares a tuning variable named `name` with uniqueId `id` and all the semantic information stored in `info`. Note that the VariableInfo struct has a `void*` field called `toolProvidedInfo`. If you fill this in, every time you get a value of that type you'll also get back that same pointer.

```c++
void kokkosp_declare_input_type(const char*, const size_t, Kokkos_Tools_VariableInfo& info);
```

This is almost exactly like declaring a tuning variable. The only difference is that in cases where the candidate values aren't known, `info.valueQuantity` will be set to `kokkos_value_unbounded`. This is fairly common, Kokkos can tell you that `kernel_name` is a string, but we can't tell you what strings a user might provide.

```c++
void kokkosp_request_values(
    const size_t contextId,
    const size_t numContextVariables, const Kokkos_Tools_VariableValue* contextVariableValues,
    const size_t numTuningVariables, Kokkos_Tools_VariableValue* tuningVariableValues);
```

Here Kokkos is requesting the values of tuning variables, and most of the meat is here. The contextId tells us the scope across which these variables were used.

The next two arguments describe the context you're tuning in. You have the number of context variables, and an array of that size containing their values. Note that the Kokkos_Tuning_VariableValue has a field called `metadata` containing all the info (type, semantics, and critically, candidates) about that variable.

The two arguments following those describe the Tuning Variables. First the number of them, then an array of that size which you can overwrite. *Overwriting those values is how you give values back to the application*

Critically, as tuningVariableValues comes preloaded with default values, if your function body is `return;` you will not crash Kokkos, only make us use our defaults. If you don't know, you are allowed to punt and let Kokkos do what it would.

```c++
void kokkosp_begin_context(size_t contextId);
```

This starts the context pointed at by contextId. If tools use measurements to drive tuning, this is where they'll do their starting measurement.

```c++
void kokkosp_end_context(const size_t contextId);
```

This simply says that the contextId in the argument is now over. If you provided tuning values associated with that context, those values can now be associated with a result.

### App Implementation

For 99% of applications, all you need to do to interact with Kokkos Tuning Tools in your code is nothing. The only exceptions are if you want the tuning to be aware of what's happening in your application (number of particles active, whether different physics are active) if
you think that might change what the Tuning decides. If you're feeling especially brave, you can also use the Tuning interface to tune parameters within your own application. For making people aware of your application context, you need to know about a few functions


```c++
size_t Kokkos::Tools::Experimental::declare_input_type(const std::string& variableName
                            VariableInfo info,
                            );
```

This function tells a tool that you have some variable they should know about when tuning. The info describes the semantics of your variable. This is discussed in great detail under "Semantics of Variables", but you need to say whether the values will be text, int, or float, whether they're categorical, ordinal,interval, or ratio data, and whether the candidate values are "unbounded" (if you don't know the full set of values), a set, or a range. This returns a `size_t` that you should store, it's how you'll later identify what values you're providing or requesting from the tool. Note that this call doesn't actually tell the tools about values, it simply tells the tool about the nature of values you'll provide later.


```c++
size_t Kokkos::Tools::Experimental::get_new_context_id();
size_t Kokkos::Tools::Experimental::get_current_context_id();
```

    In this interface,
    you will associate values with
    "contexts" in order to decide when a given declaration of a value has gone
        out of scope.The first gets you a new context
            ID if you 're starting some new set of values. If you need to recover the last context ID so you can append to that context, rather than overwriting it with a new one, you can use `get_current_context_id()`. You' ll
                use that context id to start a context in the function

```c++ void Kokkos::Tools::Experimental::begin_context(size_t context_id);
```

This tells the tool that you're beginning a region in which you'll be setting and requesting values. If the tool optimizes for time, you're telling them to start their timer.


```c++
void Kokkos::Tools::Experimental::set_input_values(size_t contextId, size_t count,
                                  VariableValue* values);
```

Here you tell tools the values for your context variables. The contextId is used to later tell when this has gone out of scope, the count is how many variables you're declaring, and the values should come from calling `Kokkos::Tools::Experimental::make_variable_value` with the appropriate variable ID and value.

```c++
void Kokkos::Tools::Experimental::end_context(size_t contextId);
```

    This tells the tool that values from this context are no longer valid,
    and that the tool should stop their timers.

        For those who want to declare and request tuning variables,
    you only need two more functions.

```c++ void Kokkos::Tools::Experimental::declare_output_type(
        const std::string&variableName VariableInfo info);
```

    This is exactly like declareContextVariable.The only difference is that
        the
            ID's this returns should be passed to request_output_values, and that the `candidates` field in the info _must_ list valid values for the tool to provide.

```c++ void Kokkos::Tools::Experimental::request_output_values(
                size_t contextId, size_t count, VariableValue* values, );
```

Here is where you request that the tool give you a set of values. You need a contextId so that the tool can know when you're done using the value and measure results. The count tells the tool how many variables it's providing values for. Values is an array of your default values for that parameter, it must not crash your program if unchanged.

### Kokkos implementation

In the past, Kokkos and Kokkos-tools didn't share source code. Except for a "SpaceHandle" struct which users manually copied to their tools, nothing from Kokkos hit the tools repo, the interface consisted entirely of basic C types. If you read the ideas section, it translates to a lot of structs and enums. Despite my best efforts to minimize them, I think we now need to share some header files with kokkos-tools. Andrew Gaspar did really excellent work making this practical, we have

1) Kokkos_Profiling_C_Interface.h , which is (shockingly) a C interface that everything in Kokkos tools boils down to
2) Kokkos_Profiling_Interface.hpp, nice C++ wrappers around the C so that the C idioms don't hit Kokkos
3) Kokkos_Profiling.[cpp/hpp], which contain things Kokkos needs to implement tooling, but the tools don't need to know about

All of our function pointer initialization and all that mess now go into Kokkos_Profiling.[cpp/hpp], all the types are in the Interface files. The interface files will be shared with kokkos/kokkos-tools.

In terms of build changes, we now install the above .h file, and have a KOKKOS_ENABLE_TUNING option
