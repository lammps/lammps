/*
//@HEADER
// ************************************************************************
//
//                        Kokkos v. 3.0
//       Copyright (2020) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY NTESS "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL NTESS OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Christian R. Trott (crtrott@sandia.gov)
//
// ************************************************************************
//@HEADER
*/
/**
 * Before digging in to the code, it's worth taking a moment to review this
 * design. Fundamentally, what we're looking to do is allow people to test that
 * a piece of code produces some expected series of tool events. Maybe we want
 * to check that deep_copy on an execution space instance only causes the
 * expected types of fences, or that calls to resize(WithoutInitializing,...)
 * don't call an initialization kernel.
 *
 * This design is realized with an interface in which you provide a code region,
 * and a set of matchers that consume the events. These matchers are lambdas
 * that accept some set of tool events, and analyze their content, and return
 * success or failure.
 *
 * Digging into implementation, this works by having a class hierarchy of Tool
 * Events, rooted at EventBase. Every Tool event inherits from this
 * (BeginParallelForEvent, PushRegionEvent, etc). We subscribe a Kokkos Tool
 * that pushes instances of these events into a vector as a code region runs. We
 * then iterate over the list of events and the matchers, first making sure that
 * every event is of the right type to be used in the matcher, and then passing
 * it to the matcher.
 *
 * Current examples are in TestEventCorrectness.hpp
 */

#include <Kokkos_Core.hpp>
#include <sstream>
#include <iostream>
#include <utility>
#include <type_traits>
namespace Kokkos {

namespace Test {

namespace Tools {

/**
 * @brief This is what a matcher should return
 * It is a two-part struct, with a bool representing
 * success (true if the match holds), and a vector of
 * strings representing the diagnostics that should be
 * printed in case of a failure
 */
struct MatchDiagnostic {
  bool success                      = true;
  std::vector<std::string> messages = {};
};

struct EventBase;  // forward declaration
using EventBasePtr = std::shared_ptr<EventBase>;
using event_vector = std::vector<EventBasePtr>;

/**
 * @brief Base case of a recursive reduction using templates
 * Should be replaced with a fold in C++17
 */

inline bool are_valid() { return true; }

/**
 * @brief Recursive reduction to check whether any pointer in a set is null
 *
 * @tparam Head Type of the pointer to examine
 * @tparam Tail Types of the rest of the pointers
 * @param head The pointer to examine
 * @param tail The rest of the pointers
 * @return true if no pointer is null, false otherwise
 *
 */
template <class Head, class... Tail>
bool are_valid(const Head& head, const Tail&... tail) {
  return (head != nullptr) && (are_valid(tail...));
}

/**
 * @brief In order to call some arbitrary set of lambdas representing matchers,
 * we need the ability to look at a lambda, and deduce its arguments.
 *
 * This is the base template, and will be specialized. All specializations
 * should define
 * - a return type R,
 * - an args pack A,
 * - a num_args, and
 * - a function "invoke_as" that takes a functor and an arg-pack, and tries to
 * call the functor with that arg-pack.
 *
 * The main original intent here is two-fold, one to allow us to look at how
 * many args a functor takes, and two to look at the types of its args. The
 * second of these is used to do a series of dynamic_casts, making sure that the
 * EventBase instances captured in our event vectors are of the types being
 * looked for by our matchers
 *
 * @tparam T a functor-like object
 * @tparam typename used for specialization shenanigans
 */
template <typename T, typename = void>
struct function_traits;

/**
 * @brief Specialization of function traits, representing a free function.
 * See the base template for info on what this struct is doing.
 *
 * @tparam R return type of the function
 * @tparam A arg pack
 */
template <typename R, typename... A>
struct function_traits<R (*)(A...)> {
  using return_type                  = R;
  using class_type                   = void;
  using args_type                    = std::tuple<A...>;
  constexpr static int num_arguments = sizeof...(A);
  template <class Call, class... Args>
  static auto invoke_as(const Call& call, Args&&... args) {
    if (!are_valid(std::dynamic_pointer_cast<A>(std::forward<Args>(args))...)) {
      return MatchDiagnostic{false, {"Types didn't match on arguments"}};
    }
    return call(*std::dynamic_pointer_cast<A>(std::forward<Args>(args))...);
  }
};

/**
 * @brief Specialization of function traits, representing a class member
 * function. See the base template for info on what this struct is doing
 *
 * @tparam R return type of the function
 * @tparam C the class function being represented
 * @tparam A arg pack
 */

template <typename R, typename C, typename... A>
struct function_traits<R (C::*)(A...)> {
  using return_type                  = R;
  using class_type                   = void;
  using args_type                    = std::tuple<A...>;
  constexpr static int num_arguments = sizeof...(A);
  template <class Call, class... Args>
  static auto invoke_as(const Call& call, Args&&... args) {
    if (!are_valid(std::dynamic_pointer_cast<A>(std::forward<Args>(args))...)) {
      return MatchDiagnostic{false, {"Types didn't match on arguments"}};
    }
    return call(*std::dynamic_pointer_cast<A>(std::forward<Args>(args))...);
  }
};

/**
 * @brief Specialization of function traits, representing a *const* class member
 * function. See the base template for info on what this struct is doing
 *
 * @tparam R return type of the function
 * @tparam C the class function being represented
 * @tparam A arg pack
 */

template <typename R, typename C, typename... A>
struct function_traits<R (C::*)(A...) const>  // const
{
  using return_type                  = R;
  using class_type                   = C;
  using args_type                    = std::tuple<A...>;
  constexpr static int num_arguments = sizeof...(A);
  template <class Call, class... Args>
  static auto invoke_as(const Call& call, Args&&... args) {
    if (!are_valid(std::dynamic_pointer_cast<A>(std::forward<Args>(args))...)) {
      return MatchDiagnostic{false, {"Types didn't match on arguments"}};
    }
    return call(*std::dynamic_pointer_cast<A>(std::forward<Args>(args))...);
  }
};

/**
 * @brief Specialization of function traits, representing a T that has a
 * non-generic call operator, i.e. a functor/lambda whose operator() has no auto
 * or template on it. See the base template for info on what this struct is
 * doing.
 *
 * @tparam T The functor type
 */
template <typename T>
struct function_traits<T, Kokkos::Impl::void_t<decltype(&T::operator())> >
    : public function_traits<decltype(&T::operator())> {};

/**
 * @brief A struct to extract events from an event vector, and invoke a matcher
 * with them.
 *
 * This one is a bit funky, you can't do std::get's or the like with a vector.
 * So this takes in a number of arguments to pull from the vector, and a start
 * index at which to begin taking from. It then makes an index sequence of that
 * number of elements {0, 1, 2, ..., num}, and then uses the function_traits
 * trick above to invoke the matcher with
 * {events[index+0],events[index+1],...,events[index+num-1]}.
 *
 * @tparam num number of arguments to the functor
 * @tparam Matcher the lambda we want to call with events from our event vector
 */
template <int num, class Matcher>
struct invoke_helper {
 private:
  // private helper with an index_sequence, invokes the matcher
  template <class Traits, size_t... Indices>
  static auto call(int index, const event_vector& events,
                   std::index_sequence<Indices...>, const Matcher& matcher) {
    return Traits::invoke_as(matcher, events[index + Indices]...);
  }

 public:
  // the entry point to the class, takes in a Traits class that knows how to
  // invoke the matcher,
  template <class Traits>
  static auto call(int index, const event_vector& events,
                   const Matcher& matcher) {
    return call<Traits>(index, events, std::make_index_sequence<num>{},
                        matcher);
  }
};

/**
 * @brief This is the base case of a recursive check of matchers, meaning no
 * more matchers exist. The only check now should be that we made it all the way
 * through the list of events captured by our lambda.
 *
 * @param events_scanned how many events we scanned
 * @param events the vector containing our events
 * @return MatchDiagnostic success if we scanned all events, failure otherwise
 */
inline MatchDiagnostic check_match(event_vector::size_type events_scanned,
                                   const event_vector& events) {
  auto result =
      ((events_scanned == events.size())
           ? MatchDiagnostic{true}
           : MatchDiagnostic{false, {"Wrong number of events encountered"}});
  return result;
}

/**
 * @brief Checks that a set of matchers match the events produced by a code
 * region
 *
 * @tparam Matcher a functor that accepts a set of events, and returns whether
 * they meet an expected structure
 * @tparam Matchers additional matchers to invoke, supposing the current one is
 * fine
 * @param index What position in our vector of events to begin pulling events
 * from
 * @param events A vector of events we want to match against our matchers
 * @param matcher the instance of Matcher (see above)
 * @param matchers the instances of Matchers (see above)
 * @return MatchDiagnostic success if the matcher matches, failure otherwise
 */
template <class Matcher, class... Matchers>
MatchDiagnostic check_match(event_vector::size_type index,
                            const event_vector& events, const Matcher& matcher,
                            const Matchers&... matchers) {
  // struct that tells us what we want to know about our matcher, and helps us
  // invoke it
  using Traits = function_traits<Matcher>;
  // how many args does that lambda have?
  constexpr static event_vector::size_type num_args = Traits::num_arguments;
  // make sure that we don't have insufficient events in our event vector
  if (index + num_args > events.size()) {
    return {false, {"Not enough events encountered to fill the matchers"}};
  }
  // Call the lambda, if it's callable with our args. Store the resulting
  // MatchDiagnostic
  auto result = invoke_helper<num_args, Matcher>::template call<Traits>(
      index, events, matcher);
  // If we fail, don't continue looking for more matches, just return
  if (!result.success) {
    return result;
  }
  // Otherwise, call with the next matcher
  return check_match(index + num_args, events, matchers...);
}

/**
 * @brief Small utility helper, an entry point into "check_match."
 * The real "check_match" needs an index at which to start checking,
 * this just tells it "hey, start at 0."
 *
 */
template <class... Matchers>
auto check_match(const event_vector& events, Matchers&&... matchers) {
  return check_match(0, events, std::forward<Matchers>(matchers)...);
}

/**
 * @brief Base class of representing everything you can do with an Event
 * checked by this system. Not much is required, just the ability to
 * represent yourself as a string for debugging purposes
 */
struct EventBase {
  using PtrHandle                        = const void* const;
  virtual ~EventBase()                   = default;
  virtual std::string descriptor() const = 0;
};

/**
 * @brief There are an unholy number of begin events in Kokkos, this is a base
 * class for them (BeginParallel[For/Reduce/Scan], BeginFence).
 *
 * @tparam Derived CRTP, intended for use with dynamic_casts
 */
template <class Derived>
struct BeginOperation : public EventBase {
  const std::string name;
  const uint32_t deviceID;
  uint64_t kID;
  BeginOperation(const std::string& n, const uint32_t devID, uint64_t k)
      : name(n), deviceID(devID), kID(k) {}
  std::string descriptor() const override {
    std::stringstream s;
    s << Derived::begin_op_name() << " { \"" << name << "\", ";
    s << deviceID;
    s << ",";
    s << kID;
    s << "}";
    return s.str();
  }
};
/**
 * @brief Analogous to BeginOperation, there are a lot of things in Kokkos
 * of roughly this structure.
 *
 * @tparam Derived CRTP, used for comparing that EventBase instances are of the
 * same type
 */
template <class Derived>
struct EndOperation : public EventBase {
  uint64_t kID;
  EndOperation(uint64_t k) : kID(k) {}

  std::string descriptor() const override {
    std::stringstream s;
    s << Derived::end_op_name() << " { ";
    s << kID;
    s << "}";
    return s.str();
  }
};

/**
 * Note, the following classes look identical, and they are. They exist because
 * we're using dynamic_casts up above to check whether events are of the same
 * type. So the different type names here are meaningful, even though the
 * classes are empty
 */
struct BeginParallelForEvent : public BeginOperation<BeginParallelForEvent> {
  static const std::string& begin_op_name() {
    static std::string value = "BeginParallelFor";
    return value;
  }
  BeginParallelForEvent(std::string n, const uint32_t devID, uint64_t k)
      : BeginOperation<BeginParallelForEvent>(n, devID, k) {}
};
struct BeginParallelReduceEvent
    : public BeginOperation<BeginParallelReduceEvent> {
  static const std::string& begin_op_name() {
    static std::string value = "BeginParallelReduce";
    return value;
  }

  BeginParallelReduceEvent(std::string n, const uint32_t devID, uint64_t k)
      : BeginOperation<BeginParallelReduceEvent>(n, devID, k) {}
};
struct BeginParallelScanEvent : public BeginOperation<BeginParallelScanEvent> {
  static const std::string& begin_op_name() {
    static std::string value = "BeginParallelScan";
    return value;
  }

  BeginParallelScanEvent(std::string n, const uint32_t devID, uint64_t k)
      : BeginOperation<BeginParallelScanEvent>(n, devID, k) {}
};
struct BeginFenceEvent : public BeginOperation<BeginFenceEvent> {
  static const std::string& begin_op_name() {
    static std::string value = "BeginFence";
    return value;
  }

  BeginFenceEvent(std::string n, const uint32_t devID, uint64_t k)
      : BeginOperation<BeginFenceEvent>(n, devID, k) {}
};

struct EndParallelForEvent : public EndOperation<EndParallelForEvent> {
  static const std::string& end_op_name() {
    static std::string value = "EndParallelFor";
    return value;
  }

  EndParallelForEvent(uint64_t k) : EndOperation<EndParallelForEvent>(k) {}
};
struct EndParallelReduceEvent : public EndOperation<EndParallelReduceEvent> {
  static const std::string& end_op_name() {
    static std::string value = "EndParallelReduce";
    return value;
  }

  EndParallelReduceEvent(uint64_t k)
      : EndOperation<EndParallelReduceEvent>(k) {}
};
struct EndParallelScanEvent : public EndOperation<EndParallelScanEvent> {
  static const std::string& end_op_name() {
    static std::string value = "EndParallelScan";
    return value;
  }

  EndParallelScanEvent(uint64_t k) : EndOperation<EndParallelScanEvent>(k) {}
};
struct EndFenceEvent : public EndOperation<EndFenceEvent> {
  static const std::string& end_op_name() {
    static std::string value = "EndFence";
    return value;
  }

  EndFenceEvent(uint64_t k) : EndOperation<EndFenceEvent>(k) {}
};

struct InitEvent : public EventBase {
  int load_sequence;
  uint64_t version_number;
  uint32_t num_device_infos;
  Kokkos::Profiling::KokkosPDeviceInfo* device_infos;
  std::string descriptor() const override {
    std::stringstream s;
    s << "InitEvent { load_sequence: " << load_sequence << ", version_number "
      << version_number << ", num_device_infos " << num_device_infos << "}";
    return s.str();
  }
  InitEvent(int l, uint64_t v_n, uint32_t n_d_i,
            Kokkos::Profiling::KokkosPDeviceInfo* d_i)
      : load_sequence(l),
        version_number(v_n),
        num_device_infos(n_d_i),
        device_infos(d_i) {}
};
struct FinalizeEvent : public EventBase {
  std::string descriptor() const override { return "FinalizeEvent{}"; }
};

struct ParseArgsEvent : public EventBase {
  int num_args;
  char** args;

  std::string descriptor() const override {
    std::stringstream s;
    s << "ParseArgsEvent { num_args : " << num_args << std::endl;
    for (int x = 0; x < num_args; ++x) {
      s << "  \"" << args[x] << "\"" << std::endl;
    }
    s << "}";
    return s.str();
  }
  ParseArgsEvent(int n_a, char** a) : num_args(n_a), args(a) {}
};
struct PrintHelpEvent : public EventBase {
  char* prog_name;
  std::string descriptor() const override {
    return "PrintHelpEvent { Program Name: \"" + std::string(prog_name) + "\"}";
  }
  PrintHelpEvent(char* p_n) : prog_name(p_n) {}
};
struct PushRegionEvent : public EventBase {
  std::string name;
  std::string descriptor() const override {
    return "PushRegionEvent { Region Name: \"" + name + "\" }";
  }
  PushRegionEvent(std::string n) : name(n) {}
};
struct PopRegionEvent : public EventBase {
  std::string descriptor() const override { return "PopRegionEvent{}"; }
};

template <class Derived>
struct DataEvent : public EventBase {
  using SpaceHandleType = Kokkos::Profiling::SpaceHandle;
  SpaceHandleType handle;
  std::string name;
  EventBase::PtrHandle ptr;
  uint64_t size;

  std::string descriptor() const override {
    std::stringstream s;
    s << Derived::event_name() << "{ In space \"" << handle.name
      << "\", name: \"" << name << "\", ptr: " << ptr << ", size: " << size
      << "}";
    return s.str();
  }
  DataEvent(SpaceHandleType h, std::string n, EventBase::PtrHandle p,
            uint64_t s)
      : handle(h), name(n), ptr(p), size(s) {}
};

struct AllocateDataEvent : public DataEvent<AllocateDataEvent> {
  static std::string event_name() { return "AllocateDataEvent"; }
  AllocateDataEvent(DataEvent::SpaceHandleType h, std::string n,
                    EventBase::PtrHandle p, uint64_t s)
      : DataEvent<AllocateDataEvent>(h, n, p, s) {}
};
struct DeallocateDataEvent : public DataEvent<DeallocateDataEvent> {
  static std::string event_name() { return "DeallocateDataEvent"; }
  DeallocateDataEvent(DataEvent::SpaceHandleType h, std::string n,
                      EventBase::PtrHandle p, uint64_t s)
      : DataEvent<DeallocateDataEvent>(h, n, p, s) {}
};

struct CreateProfileSectionEvent : public EventBase {
  std::string name;
  uint32_t id;
  std::string descriptor() const override {
    return "CreateProfileSectionEvent {\"" + name + "\", " +
           std::to_string(id) + "}";
  }
  CreateProfileSectionEvent(std::string n, uint32_t s_i) : name(n), id(s_i) {}
};

template <class Derived>
struct ProfileSectionManipulationEvent : public EventBase {
  uint32_t id;
  std::string descriptor() const override {
    std::stringstream s;
    s << Derived::event_name() << "{ " << id << "}";
    return s.str();
  }
  ProfileSectionManipulationEvent(uint32_t d_i) : id(d_i){};
};

struct StartProfileSectionEvent
    : public ProfileSectionManipulationEvent<StartProfileSectionEvent> {
  static std::string event_name() { return "StartProfileSectionEvent"; }
  StartProfileSectionEvent(uint32_t d_i)
      : ProfileSectionManipulationEvent<StartProfileSectionEvent>(d_i){};
};
struct StopProfileSectionEvent
    : public ProfileSectionManipulationEvent<StopProfileSectionEvent> {
  static std::string event_name() { return "StopProfileSectionEvent"; }
  StopProfileSectionEvent(uint32_t d_i)
      : ProfileSectionManipulationEvent<StopProfileSectionEvent>(d_i){};
};
struct DestroyProfileSectionEvent
    : public ProfileSectionManipulationEvent<DestroyProfileSectionEvent> {
  static std::string event_name() { return "DestroyProfileSectionEvent"; }
  DestroyProfileSectionEvent(uint32_t d_i)
      : ProfileSectionManipulationEvent<DestroyProfileSectionEvent>(d_i){};
};

struct ProfileEvent : public EventBase {
  std::string name;
  std::string descriptor() const override {
    return "ProfileEvent {\"" + name + "\"}";
  }
  ProfileEvent(std::string n) : name(n) {}
};

struct BeginDeepCopyEvent : public EventBase {
  using SpaceHandleType = Kokkos::Profiling::SpaceHandle;
  SpaceHandleType src_handle;
  std::string src_name;
  EventBase::PtrHandle src_ptr;
  SpaceHandleType dst_handle;
  std::string dst_name;
  EventBase::PtrHandle dst_ptr;
  uint64_t size;
  std::string descriptor() const override {
    std::stringstream s;
    s << "BeginDeepCopyEvent { size: " << size << std::endl;
    s << "  dst: { \"" << dst_handle.name << "\", \"" << dst_name << "\", "
      << dst_ptr << "}\n";
    s << "  src: { \"" << src_handle.name << "\", \"" << src_name << "\", "
      << src_ptr << "}\n";
    s << "}";
    return s.str();
  }
  BeginDeepCopyEvent(SpaceHandleType s_h, std::string s_n,
                     EventBase::PtrHandle s_p, SpaceHandleType d_h,
                     std::string d_n, EventBase::PtrHandle d_p, uint64_t s)
      : src_handle(s_h),
        src_name(s_n),
        src_ptr(s_p),
        dst_handle(d_h),
        dst_name(d_n),
        dst_ptr(d_p),
        size(s) {}
};
struct EndDeepCopyEvent : public EventBase {
  std::string descriptor() const override { return "EndDeepCopyEvent{}"; }
};

template <class Derived>
struct DualViewEvent : public EventBase {
  std::string name;
  EventBase::PtrHandle ptr;
  bool is_device;
  DualViewEvent(std::string n, EventBase::PtrHandle p, bool i_d)
      : name(n), ptr(p), is_device(i_d) {}
  std::string descriptor() const override {
    std::stringstream s;
    s << Derived::event_name() << " { \"" << name << "\", " << std::hex << ptr
      << ", " << std::boolalpha << is_device << "}";
    return s.str();
  }
};
struct DualViewModifyEvent : public DualViewEvent<DualViewModifyEvent> {
  static std::string event_name() { return "DualViewModifyEvent"; }
  DualViewModifyEvent(std::string n, EventBase::PtrHandle p, bool i_d)
      : DualViewEvent(n, p, i_d) {}
};
struct DualViewSyncEvent : public DualViewEvent<DualViewSyncEvent> {
  static std::string event_name() { return "DualViewSyncEvent"; }
  DualViewSyncEvent(std::string n, EventBase::PtrHandle p, bool i_d)
      : DualViewEvent(n, p, i_d) {}
};

struct DeclareMetadataEvent : public EventBase {
  std::string key;
  std::string value;
  std::string descriptor() const override {
    return "DeclareMetadataEvent {\"" + key + "\", \"" + value + "\"}";
  }
  DeclareMetadataEvent(std::string k, std::string v) : key(k), value(v) {}
};

struct ProvideToolProgrammingInterfaceEvent : public EventBase {
  using Interface = Kokkos::Tools::Experimental::ToolProgrammingInterface;

  uint32_t num_functions;
  Interface interface;
  ProvideToolProgrammingInterfaceEvent(uint32_t n_f, Interface i)
      : num_functions(n_f), interface(i) {}
  std::string descriptor() const override {
    return "ProvideToolProgrammingInterfaceEvent {" +
           std::to_string(num_functions) + "}";
  }
};
struct RequestToolSettingsEvent : public EventBase {
  using Settings = Kokkos::Tools::Experimental::ToolSettings;

  uint32_t num_settings;
  Settings settings;
  RequestToolSettingsEvent(uint32_t n_s, Settings s)
      : num_settings(n_s), settings(s) {}
  std::string descriptor() const override {
    return "RequestToolSettingsEvent {" + std::to_string(num_settings) + "}";
  }
};

template <class Derived>
struct TypeDeclarationEvent : public EventBase {
  std::string name;
  size_t variable_id;
  Kokkos::Tools::Experimental::VariableInfo info;
  std::string descriptor() const override {
    return Derived::event_name() + "{ \"" + name + "\"," +
           std::to_string(variable_id) + "}";
  }
  TypeDeclarationEvent(std::string n, size_t v_i,
                       Kokkos::Tools::Experimental::VariableInfo i)
      : name(n), variable_id(v_i), info(i) {}
};
struct DeclareOutputTypeEvent
    : public TypeDeclarationEvent<DeclareOutputTypeEvent> {
  static std::string event_name() { return "DeclarateOutputTypeEvent"; }
  DeclareOutputTypeEvent(std::string n, size_t v_i,
                         Kokkos::Tools::Experimental::VariableInfo i)
      : TypeDeclarationEvent(n, v_i, i) {}
};
struct DeclareInputTypeEvent
    : public TypeDeclarationEvent<DeclareInputTypeEvent> {
  static std::string event_name() { return "DeclareInputTypeEvent"; }
  DeclareInputTypeEvent(std::string n, size_t v_i,
                        Kokkos::Tools::Experimental::VariableInfo i)
      : TypeDeclarationEvent(n, v_i, i) {}
};

struct RequestOutputValuesEvent : public EventBase {
  size_t context;
  size_t num_inputs;
  std::vector<Kokkos::Tools::Experimental::VariableValue> inputs;
  size_t num_outputs;
  std::vector<Kokkos::Tools::Experimental::VariableValue> outputs;
  std::string descriptor() const override {
    std::stringstream s;
    s << "RequestOutputValuesEvent { ";
    s << num_inputs << " inputs,";
    s << num_outputs << " outputs}";
    return s.str();
  }
  RequestOutputValuesEvent(
      size_t c, size_t n_i,
      std::vector<Kokkos::Tools::Experimental::VariableValue> i, size_t n_o,
      std::vector<Kokkos::Tools::Experimental::VariableValue> o)
      : context(c), num_inputs(n_i), inputs(i), num_outputs(n_o), outputs(o) {}
};

struct BeginContextEvent : public EventBase {
  size_t context;
  std::string descriptor() const override {
    return "ContextBeginEvent{ " + std::to_string(context) + "}";
  }
  BeginContextEvent(size_t c) : context(c) {}
};
struct EndContextEvent : public EventBase {
  size_t context;
  Kokkos::Tools::Experimental::VariableValue value;
  std::string descriptor() const override {
    return "ContextEndEvent {" + std::to_string(context) + "}";
  }
  EndContextEvent(size_t c, Kokkos::Tools::Experimental::VariableValue v)
      : context(c), value(v) {}
};

struct OptimizationGoalDeclarationEvent : public EventBase {
  size_t context;
  Kokkos::Tools::Experimental::OptimizationGoal goal;
  std::string descriptor() const override {
    return "OptimizationGoalDeclarationEvent{" + std::to_string(context) + "}";
  }
  OptimizationGoalDeclarationEvent(
      size_t c, Kokkos::Tools::Experimental::OptimizationGoal g)
      : context(c), goal(g) {}
};

/**
 * @brief Takes a vector of events, a set of matchers, and checks whether
 *        that event vector matches what those matchers expect
 *
 * @tparam Matchers types of our matchers
 * @param events A vector containing events
 * @param matchers A set of functors that match those Events
 * @return true on successful match, false otherwise
 */
template <class... Matchers>
bool compare_event_vectors(const event_vector& events, Matchers&&... matchers) {
  // leans on check_match to do the bulk of the work
  auto diagnostic = check_match(events, std::forward<Matchers>(matchers)...);
  // On failure, print out the error messages
  if (!diagnostic.success) {
    for (const auto& message : diagnostic.messages) {
      std::cerr << "Error matching event vectors: " << message << std::endl;
    }
  }
  return diagnostic.success;
}

/**
 * This section is odd, and needs explanation. Imagine that
 * you're writing a test. Maybe you want to listen to all
 * events. Maybe you want to listen to all profiling events.
 * Maybe you want to listen to all profiling events, no
 * infrastructure events, and only type declaration events
 * in tuning.
 *
 * You can model this as a tree of preferences, a kind of
 * hierarchical bool. By default,
 * we listen to everything. But you can disable everything,
 * or any subcomponent (profiling/tuning/infrastructure),
 * or even a sub-subcomponent (profiling->kernels)
 *
 */

/**
 * @brief This tells the testing tool which events to listen to.
 * My strong recommendation is to make this "all events" in most cases,
 * but if there is an event that is hard to match in some cases, a stray
 * deep_copy or the like, this will let you ignore that event. Users will
 * not directly instantiate these.
 */

struct ToolValidatorConfiguration {
  struct Profiling {
    bool kernels        = true;
    bool regions        = true;
    bool fences         = true;
    bool allocs         = true;
    bool copies         = true;
    bool dual_view_ops  = true;
    bool sections       = true;
    bool profile_events = true;
    bool metadata       = true;
  };
  struct Tuning {
    bool contexts          = true;
    bool type_declarations = true;
    bool request_values    = true;
  };
  struct Infrastructure {
    bool init                  = true;
    bool finalize              = true;
    bool programming_interface = true;
    bool request_settings      = true;
  };
  Profiling profiling           = Profiling();
  Tuning tuning                 = Tuning();
  Infrastructure infrastructure = Infrastructure();
};

namespace Config {
/**
 * @brief A config struct has a few properties:
 *
 * 1) What settings it toggles
 * 2) Whether it toggles that setting on or off
 * 3) What depth the setting is in the tree
 *
 * The first two hopefully make intuitive sense. The
 * third is weird. In order to make this hierarchical
 * bool concept work, you need to be able to first
 * disable all events, then enable profiling.
 *
 * This is done by modeling the depth of the request.
 * DisableAlls happen before EnableProfiling happen before
 * DisableKernels. The implementation of that is in listen_tool_events,
 * but needs machinery here.
 *
 */

/**
 * @brief Macro to make defining a configuration struct easier.
 * Given a name, what value to override in the ToolConfiguration,
 * and the depth of that configuration option, produces an
 * EnableName struct to enable that option, and a DisableName
 * struct to disable that option
 *
 * @param name : the name of the struct
 * @param value: the value in ToolConfiguration to override
 * @param depth: how deep in the configuration tree an option is
 *               (0 is root, Profiling/Tuning/Infrastructure 1, 2 for
 *                sub-options)
 */
#define KOKKOS_IMPL_TOOLS_TEST_CONFIG_OPTION(name, value, depth)    \
  template <bool target_value>                                      \
  struct Toggle##name : public std::integral_constant<int, depth> { \
    void operator()(ToolValidatorConfiguration& config) {           \
      config.value = target_value;                                  \
    }                                                               \
  };                                                                \
  using Enable##name  = Toggle##name<true>;                         \
  using Disable##name = Toggle##name<false>

KOKKOS_IMPL_TOOLS_TEST_CONFIG_OPTION(Kernels, profiling.kernels, 2);
KOKKOS_IMPL_TOOLS_TEST_CONFIG_OPTION(Regions, profiling.regions, 2);
KOKKOS_IMPL_TOOLS_TEST_CONFIG_OPTION(Fences, profiling.fences, 2);
KOKKOS_IMPL_TOOLS_TEST_CONFIG_OPTION(Allocs, profiling.allocs, 2);
KOKKOS_IMPL_TOOLS_TEST_CONFIG_OPTION(Copies, profiling.copies, 2);
KOKKOS_IMPL_TOOLS_TEST_CONFIG_OPTION(DualViewOps, profiling.dual_view_ops, 2);
KOKKOS_IMPL_TOOLS_TEST_CONFIG_OPTION(Sections, profiling.sections, 2);
KOKKOS_IMPL_TOOLS_TEST_CONFIG_OPTION(ProfileEvents, profiling.profile_events,
                                     2);
KOKKOS_IMPL_TOOLS_TEST_CONFIG_OPTION(Metadata, profiling.metadata, 2);
KOKKOS_IMPL_TOOLS_TEST_CONFIG_OPTION(Contexts, tuning.contexts, 2);
KOKKOS_IMPL_TOOLS_TEST_CONFIG_OPTION(TypeDeclarations, tuning.type_declarations,
                                     2);
KOKKOS_IMPL_TOOLS_TEST_CONFIG_OPTION(RequestValues, tuning.request_values, 2);
KOKKOS_IMPL_TOOLS_TEST_CONFIG_OPTION(Init, infrastructure.init, 2);
KOKKOS_IMPL_TOOLS_TEST_CONFIG_OPTION(Finalize, infrastructure.finalize, 2);
KOKKOS_IMPL_TOOLS_TEST_CONFIG_OPTION(ProgrammingInterface,
                                     infrastructure.programming_interface, 2);
KOKKOS_IMPL_TOOLS_TEST_CONFIG_OPTION(RequestSettings,
                                     infrastructure.request_settings, 2);

template <bool target_value>
struct ToggleInfrastructure : public std::integral_constant<int, 1> {
  void operator()(ToolValidatorConfiguration& config) {
    ToggleInit<target_value>{}(config);
    ToggleFinalize<target_value>{}(config);
    ToggleProgrammingInterface<target_value>{}(config);
    ToggleRequestSettings<target_value>{}(config);
  }
};

using EnableInfrastructure  = ToggleInfrastructure<true>;
using DisableInfrastructure = ToggleInfrastructure<false>;

template <bool target_value>
struct ToggleProfiling : public std::integral_constant<int, 1> {
  void operator()(ToolValidatorConfiguration& config) {
    ToggleKernels<target_value>{}(config);
    ToggleRegions<target_value>{}(config);
    ToggleFences<target_value>{}(config);
    ToggleAllocs<target_value>{}(config);
    ToggleCopies<target_value>{}(config);
    ToggleDualViewOps<target_value>{}(config);
    ToggleSections<target_value>{}(config);
    ToggleProfileEvents<target_value>{}(config);
    ToggleMetadata<target_value>{}(config);
  }
};

using EnableProfiling  = ToggleProfiling<true>;
using DisableProfiling = ToggleProfiling<false>;

template <bool target_value>
struct ToggleTuning : public std::integral_constant<int, 1> {
  void operator()(ToolValidatorConfiguration& config) {
    ToggleContexts<target_value>{}(config);
    ToggleTypeDeclarations<target_value>{}(config);
    ToggleRequestValues<target_value>{}(config);
  }
};

using EnableTuning  = ToggleTuning<true>;
using DisableTuning = ToggleTuning<false>;

template <bool target_value>
struct ToggleAll : public std::integral_constant<int, 0> {
  void operator()(ToolValidatorConfiguration& config) {
    ToggleProfiling<target_value>{}(config);
    ToggleTuning<target_value>{}(config);
    ToggleInfrastructure<target_value>{}(config);
  }
};

using EnableAll  = ToggleAll<true>;
using DisableAll = ToggleAll<false>;
}  // namespace Config

/**
 * This is the vector tool callbacks will push events into.
 * It needs to be outside of functions (to be global) because
 * it needs to be used in the tools callbacks, which are function pointers,
 * which can't capture variables. Thus we need something that doesn't require
 * capturing. In short, a global variable. :(
 */
static std::vector<EventBasePtr> found_events;
/**
 * Needs to stand outside of functions, this is the kID of the last encountered
 * begin event
 */
static uint64_t last_kernel_id;
/**
 * Needs to stand outside of functions, this is the section ID of the last
 * encountered section id
 */
static uint32_t last_section_id;

/** Subscribes to all of the requested callbacks */
static void set_tool_events_impl(const ToolValidatorConfiguration& config) {
  Kokkos::Tools::Experimental::pause_tools();  // remove all events
  if (config.profiling.kernels) {
    Kokkos::Tools::Experimental::set_begin_parallel_for_callback(
        [](const char* n, const uint32_t d, uint64_t* k) {
          *k = ++last_kernel_id;
          found_events.push_back(
              std::make_shared<BeginParallelForEvent>(std::string(n), d, *k));
        });
    Kokkos::Tools::Experimental::set_begin_parallel_reduce_callback(
        [](const char* n, const uint32_t d, uint64_t* k) {
          *k = ++last_kernel_id;
          found_events.push_back(std::make_shared<BeginParallelReduceEvent>(
              std::string(n), d, *k));
        });
    Kokkos::Tools::Experimental::set_begin_parallel_scan_callback(
        [](const char* n, const uint32_t d, uint64_t* k) {
          *k = ++last_kernel_id;

          found_events.push_back(
              std::make_shared<BeginParallelScanEvent>(std::string(n), d, *k));
        });
    Kokkos::Tools::Experimental::set_end_parallel_for_callback(
        [](const uint64_t k) {
          found_events.push_back(std::make_shared<EndParallelForEvent>(k));
        });
    Kokkos::Tools::Experimental::set_end_parallel_reduce_callback(
        [](const uint64_t k) {
          found_events.push_back(std::make_shared<EndParallelReduceEvent>(k));
        });
    Kokkos::Tools::Experimental::set_end_parallel_scan_callback(
        [](const uint64_t k) {
          found_events.push_back(std::make_shared<EndParallelScanEvent>(k));
        });
  }  // if profiling.kernels
  if (config.profiling.regions) {
    Kokkos::Tools::Experimental::set_push_region_callback([](const char* name) {
      found_events.push_back(
          std::make_shared<PushRegionEvent>(std::string(name)));
    });
    Kokkos::Tools::Experimental::set_pop_region_callback(
        []() { found_events.push_back(std::make_shared<PopRegionEvent>()); });
  }
  if (config.profiling.fences) {
    Kokkos::Tools::Experimental::set_begin_fence_callback(
        [](const char* n, const uint32_t d, uint64_t* k) {
          *k = ++last_kernel_id;
          found_events.push_back(
              std::make_shared<BeginFenceEvent>(std::string(n), d, *k));
        });

    Kokkos::Tools::Experimental::set_end_fence_callback([](const uint64_t k) {
      found_events.push_back(std::make_shared<EndFenceEvent>(k));
    });
  }  // profiling.fences
  if (config.profiling.allocs) {
    Kokkos::Tools::Experimental::set_allocate_data_callback(
        [](Kokkos::Tools::SpaceHandle handle, const char* name,
           EventBase::PtrHandle const ptr, const uint64_t size) {
          found_events.push_back(std::make_shared<AllocateDataEvent>(
              handle, std::string(name), ptr, size));
        });
    Kokkos::Tools::Experimental::set_deallocate_data_callback(
        [](Kokkos::Tools::SpaceHandle handle, const char* name,
           EventBase::PtrHandle const ptr, const uint64_t size) {
          found_events.push_back(std::make_shared<DeallocateDataEvent>(
              handle, std::string(name), ptr, size));
        });
  }
  if (config.profiling.copies) {
    Kokkos::Tools::Experimental::set_begin_deep_copy_callback(
        [](Kokkos::Tools::SpaceHandle dst_handle, const char* dst_name,
           EventBase::PtrHandle dst_ptr, Kokkos::Tools::SpaceHandle src_handle,
           const char* src_name, EventBase::PtrHandle src_ptr, uint64_t size) {
          found_events.push_back(std::make_shared<BeginDeepCopyEvent>(
              dst_handle, std::string(dst_name), dst_ptr, src_handle,
              std::string(src_name), src_ptr, size));
        });
    Kokkos::Tools::Experimental::set_end_deep_copy_callback(
        []() { found_events.push_back(std::make_shared<EndDeepCopyEvent>()); });
  }
  if (config.profiling.dual_view_ops) {
    Kokkos::Tools::Experimental::set_dual_view_sync_callback(
        [](const char* name, EventBase::PtrHandle ptr, bool is_device) {
          found_events.push_back(std::make_shared<DualViewSyncEvent>(
              std::string(name), ptr, is_device));
        });
    Kokkos::Tools::Experimental::set_dual_view_modify_callback(
        [](const char* name, EventBase::PtrHandle ptr, bool is_device) {
          found_events.push_back(std::make_shared<DualViewModifyEvent>(
              std::string(name), ptr, is_device));
        });
  }
  if (config.profiling.sections) {
    Kokkos::Tools::Experimental::set_create_profile_section_callback(
        [](const char* name, uint32_t* id) {
          *id = (++last_section_id);
          found_events.push_back(std::make_shared<CreateProfileSectionEvent>(
              std::string(name), *id));
        });
    Kokkos::Tools::Experimental::set_destroy_profile_section_callback(
        [](uint32_t id) {
          found_events.push_back(
              std::make_shared<DestroyProfileSectionEvent>(id));
        });
    Kokkos::Tools::Experimental::set_start_profile_section_callback(
        [](uint32_t id) {
          found_events.push_back(
              std::make_shared<StartProfileSectionEvent>(id));
        });
    Kokkos::Tools::Experimental::set_stop_profile_section_callback(
        [](uint32_t id) {
          found_events.push_back(std::make_shared<StopProfileSectionEvent>(id));
        });
  }
  if (config.profiling.profile_events) {
    Kokkos::Tools::Experimental::set_profile_event_callback(
        [](const char* name) {
          found_events.push_back(
              std::make_shared<ProfileEvent>(std::string(name)));
        });
  }
  if (config.profiling.metadata) {
    Kokkos::Tools::Experimental::set_declare_metadata_callback(
        [](const char* key, const char* value) {
          found_events.push_back(std::make_shared<DeclareMetadataEvent>(
              std::string(key), std::string(value)));
        });
  }
  if (config.tuning.contexts) {
    Kokkos::Tools::Experimental::set_begin_context_callback(
        [](const size_t context) {
          found_events.push_back(std::make_shared<BeginContextEvent>(context));
        });
    Kokkos::Tools::Experimental::set_end_context_callback(
        [](const size_t context,
           Kokkos::Tools::Experimental::VariableValue value) {
          found_events.push_back(
              std::make_shared<EndContextEvent>(context, value));
        });
  }
  if (config.tuning.type_declarations) {
    Kokkos::Tools::Experimental::set_declare_input_type_callback(
        [](const char* name, const size_t id,
           Kokkos::Tools::Experimental::VariableInfo* info) {
          found_events.push_back(std::make_shared<DeclareInputTypeEvent>(
              std::string(name), id, *info));
        });
    Kokkos::Tools::Experimental::set_declare_output_type_callback(
        [](const char* name, const size_t id,
           Kokkos::Tools::Experimental::VariableInfo* info) {
          found_events.push_back(std::make_shared<DeclareOutputTypeEvent>(
              std::string(name), id, *info));
        });
  }
  if (config.tuning.request_values) {
    Kokkos::Tools::Experimental::set_request_output_values_callback(
        [](const size_t context, const size_t num_inputs,
           const Kokkos::Tools::Experimental::VariableValue* inputs_in,
           const size_t num_outputs,
           Kokkos::Tools::Experimental::VariableValue* outputs_in) {
          std::vector<Kokkos::Tools::Experimental::VariableValue> inputs,
              outputs;
          std::copy(inputs_in, inputs_in + num_inputs,
                    std::back_inserter(inputs));
          std::copy(outputs_in, outputs_in + num_inputs,
                    std::back_inserter(outputs));

          found_events.push_back(std::make_shared<RequestOutputValuesEvent>(
              context, num_inputs, inputs, num_outputs, outputs));
        });
  }
  if (config.infrastructure.init) {
    Kokkos::Tools::Experimental::set_init_callback(
        [](const int loadseq, const uint64_t version, const uint32_t num_infos,
           Kokkos::Profiling::KokkosPDeviceInfo* infos) {
          found_events.push_back(
              std::make_shared<InitEvent>(loadseq, version, num_infos, infos));
        });
  }
  if (config.infrastructure.finalize) {
    Kokkos::Tools::Experimental::set_finalize_callback(
        []() { found_events.push_back(std::make_shared<FinalizeEvent>()); });
  }
  if (config.infrastructure.programming_interface) {
    Kokkos::Tools::Experimental::
        set_provide_tool_programming_interface_callback(
            [](const uint32_t num_functions,
               Kokkos::Tools::Experimental::ToolProgrammingInterface
                   interface) {
              found_events.push_back(
                  std::make_shared<ProvideToolProgrammingInterfaceEvent>(
                      num_functions, interface));
            });
  }
  if (config.infrastructure.request_settings) {
    Kokkos::Tools::Experimental::set_request_tool_settings_callback(
        [](const uint32_t num_settings,
           Kokkos::Tools::Experimental::ToolSettings* settings) {
          found_events.push_back(std::make_shared<RequestToolSettingsEvent>(
              num_settings, *settings));
        });
  }
}
template <int priority>
void listen_tool_events_impl(std::integral_constant<int, priority>,
                             ToolValidatorConfiguration&) {}

template <class Config>
void invoke_config(ToolValidatorConfiguration& in, Config conf,
                   std::true_type) {
  conf(in);
}
template <class Config>
void invoke_config(ToolValidatorConfiguration&, Config, std::false_type) {}

template <int priority, class Config, class... Configs>
void listen_tool_events_impl(std::integral_constant<int, priority> prio,
                             ToolValidatorConfiguration& in, Config conf,
                             Configs... configs) {
  invoke_config(in, conf,
                std::integral_constant<bool, priority == conf.value>{});
  listen_tool_events_impl(prio, in, configs...);
}
template <class... Configs>
void listen_tool_events(Configs... confs) {
  ToolValidatorConfiguration conf;
  listen_tool_events_impl(std::integral_constant<int, 0>{}, conf, confs...);
  listen_tool_events_impl(std::integral_constant<int, 1>{}, conf, confs...);
  listen_tool_events_impl(std::integral_constant<int, 2>{}, conf, confs...);
  set_tool_events_impl(conf);
}

/**
 * @brief This is the main entry point people will use to test their programs
 *        Given a lambda representing a code region, and a set of matchers on
 * tools events, verify that the given lambda produces events that match those
 * expected by the matchers
 *
 * @tparam Lambda Type of lam
 * @tparam Matchers Type of matchers
 * @param lam The code region that will produce events
 * @param matchers Matchers for those events, lambdas that expect events and
 * compare them
 * @return true if all events are consumed, all matchers are invoked, and all
 * matchers success, false otherwise
 */
template <class Lambda, class... Matchers>
bool validate_event_set(const Lambda& lam, Matchers&&... matchers) {
  // First, erase events from previous invocations
  found_events.clear();
  // Invoke the lambda (this will populate found_events, via tooling)
  lam();
  // compare the found events against the expected ones
  auto success =
      compare_event_vectors(found_events, std::forward<Matchers>(matchers)...);
  if (!success) {
    // on failure, print out the events we found
    for (const auto& event : found_events) {
      std::cout << event->descriptor() << std::endl;
    }
  }
  return success;
}
/**
 * @brief Analogous to validate_event_set up above, except rather than
 *        comparing to matchers, this just returns the found event vector
 *
 * @tparam Lambda as in validate_event_set
 * @param lam as in validate_event_set
 * @return auto
 */
template <class Lambda>
auto get_event_set(const Lambda& lam) {
  found_events.clear();
  lam();
  // return compare_event_vectors(expected, found_events);
  std::vector<EventBasePtr> events;
  std::copy(found_events.begin(), found_events.end(),
            std::back_inserter(events));
  return events;
}

inline MatchDiagnostic check_presence_of(const EventBasePtr&) {
  return {false};
}
template <class Matcher, class... Matchers>
MatchDiagnostic check_presence_of(const EventBasePtr& event, const Matcher& m,
                                  Matchers&&... args) {
  auto tail  = check_presence_of(event, args...);
  auto match = function_traits<Matcher>::invoke_as(m, event);
  if (tail.success) {
    for (const auto& entry : tail.messages) {
      match.messages.push_back(entry);
    }
  }
  match.success |= tail.success;
  return match;
}

template <class Lambda, class... Matchers>
bool validate_absence(const Lambda& lam, const Matchers... matchers) {
  // First, erase events from previous invocations
  found_events.clear();
  // Invoke the lambda (this will populate found_events, via tooling)
  lam();
  // compare the found events against the expected ones
  for (const auto& event : found_events) {
    MatchDiagnostic match = check_presence_of(event, matchers...);

    if (match.success) {
      std::cout << "Test failure: encountered unwanted events" << std::endl;
      for (const auto& message : match.messages) {
        std::cout << "  " << message << std::endl;
      }
      // on success, print out the events we found
      for (const auto& p_event : found_events) {
        std::cout << p_event->descriptor() << std::endl;
      }
      return false;
    }
  }
  return true;
}

}  // namespace Tools
}  // namespace Test
}  // namespace Kokkos
