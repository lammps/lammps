# This script automatically generate C++ header file
# for template multiarrays from 2 to DIMMAX dimensions

from jinja2 import Template

DIMMAX = 7

cpp_class_code_template = """template <class T> class Array{{N}}D: public ContigousArrayND<T>  {
    using ContigousArrayND<T>::array_name;
    using ContigousArrayND<T>::data;
    using ContigousArrayND<T>::size;
    
    // dimensions 
    size_t dim[{{N}}];
    
    // strides
    size_t s[{{N}}];
    
public:
    // default empty constructor
    Array{{N}}D() = default;    

    // parametrized constructor
    Array{{N}}D(const string& array_name) {this->array_name = array_name;}
   
    // parametrized constructor
    Array{{N}}D({% for i in range(N) %}size_t d{{i}}{% if not loop.last %}, {% endif %}{%endfor%}, const string& array_name="Array{{N}}D") {
        init({% for i in range(N) %}d{{i}}{% if not loop.last %}, {% endif %}{%endfor%},array_name);
    }   
    
    //initialize array and strides
    void init({% for i in range(N) %}size_t d{{i}}{% if not loop.last %}, {% endif %}{%endfor%}, const string& array_name="Array{{N}}D") {
        this->array_name = array_name;
        
        {% for i in range(N) %}
        dim[{{i}}] = d{{i}};{% endfor %}
        
        s[{{N-1}}] = 1;{% for i in range(N-2,-1,-1) %}
        s[{{i}}] = s[{{i+1}}]*dim[{{i+1}}];{% endfor %}
                
        if(size!=s[0]*dim[0]) {
            size = s[0]*dim[0];
            if(data) delete[] data;
            data = new T[size];
            memset(data,  0, size * sizeof(T));
        } else {
            memset(data,  0, size * sizeof(T));
        }
    }
    
    void resize({% for i in range(N) %}size_t d{{i}}{% if not loop.last %}, {% endif %}{%endfor%}) {
        init({% for i in range(N) %}d{{i}}{% if not loop.last %}, {% endif %}{%endfor%}, this->array_name);    
    }
    
    size_t get_dim(int d) const {
        return dim[d];
    }
  
#ifdef MULTIARRAY_INDICES_CHECK
    void check_indices({% for i in range(N) %}size_t i{{i}}{% if not loop.last %}, {% endif %}{%endfor%}) const {
    {% for i in range(N) %}
        if((i{{i}}<0)|(i{{i}}>=dim[{{i}}])){
            printf("%s: index i{{i}}=%ld out of range (0, %ld)\\n", array_name.c_str(), i{{i}}, dim[{{i}}]-1);            
            exit(EXIT_FAILURE);
        }
    {%endfor%}
    }
#endif

    inline const T& operator()({% for i in range(N) %}size_t i{{i}}{% if not loop.last %}, {% endif %}{%endfor%}) const {
#ifdef MULTIARRAY_INDICES_CHECK
        check_indices({% for i in range(N) %}i{{i}}{% if not loop.last %}, {% endif %}{%endfor%});
#endif    
{% if N > 1 %}
       return data[{% for i in range(N-1) %}i{{i}}*s[{{i}}]{% if not loop.last %}+{% endif %}{%endfor%} + i{{N-1}}];
{% else %}       
        return data[i{{N-1}}];
{% endif %}       
     }
     
    inline T& operator()({% for i in range(N) %}size_t i{{i}}{% if not loop.last %}, {% endif %}{%endfor%}) {
#ifdef MULTIARRAY_INDICES_CHECK
        check_indices({% for i in range(N) %}i{{i}}{% if not loop.last %}, {% endif %}{%endfor%});
#endif 
{% if N > 1 %}
       return data[{% for i in range(N-1) %}i{{i}}*s[{{i}}]{% if not loop.last %}+{% endif %}{%endfor%} + i{{N-1}}];
{% else %}       
        return data[i{{N-1}}];
{% endif %}           
     }
};
"""

class_template = Template(cpp_class_code_template)

cpp_general_class_code_template = """template <class T> class Array{{N}}DGeneral: public ContigousArrayND<T>  {
    using ContigousArrayND<T>::array_name;
    using ContigousArrayND<T>::data;
    using ContigousArrayND<T>::size;
    
     // ranges 
    int {% for i in range(N) %}i{{i}}_start, i{{i}}_stop {% if not loop.last %}, {% endif %}{%endfor%};
        
    // dimensions 
    size_t dim[{{N}}];
    
    // strides
    size_t s[{{N}}];
    
public:
    // default empty constructor
    Array{{N}}DGeneral() {}      
    
    // parametrized constructor
    Array{{N}}DGeneral(const string& array_name) {this->array_name = array_name;}
    
    // parametrized constructor
    Array{{N}}DGeneral({% for i in range(N) %}int i{{i}}_init, int i{{i}}_final {% if not loop.last %}, {% endif %}{%endfor%}, const string& array_name="Array{{N}}DGeneral") {
        init({% for i in range(N) %}i{{i}}_init, i{{i}}_final {% if not loop.last %}, {% endif %}{%endfor%}, array_name);
    }
    
    //initialize array and strides
    void init({% for i in range(N) %}int i{{i}}_init, int i{{i}}_final {% if not loop.last %}, {% endif %}{%endfor%}, const string& array_name="Array{{N}}DGeneral") {
        this->array_name = array_name;
        {% for i in range(N) %}
        i{{i}}_start = i{{i}}_init; i{{i}}_stop = i{{i}}_final;{% endfor %}
        
        {% for i in range(N) %}
        dim[{{i}}] = i{{i}}_final-i{{i}}_init+1;{% endfor %}
        
        s[{{N-1}}] = 1;{% for i in range(N-2,-1,-1) %}
        s[{{i}}] = s[{{i+1}}]*dim[{{i+1}}];{% endfor %}
        
        if(size!=s[0]*dim[0]) {
            size = s[0]*dim[0];
            if(data) delete[] data;
            data = new T[size];
            memset(data,  0, size * sizeof(T));
        } else {
            memset(data,  0, size * sizeof(T));
        }
    }

    void resize({% for i in range(N) %}int i{{i}}_init, int i{{i}}_final {% if not loop.last %}, {% endif %}{%endfor%}) {
        init({% for i in range(N) %}i{{i}}_init, i{{i}}_final {% if not loop.last %}, {% endif %}{%endfor%}, this->array_name);    
    }
    
    size_t get_dim(int d) const {
        return dim[d];
    }
    
#ifdef MULTIARRAY_INDICES_CHECK
    void check_indices({% for i in range(N) %}int i{{i}}{% if not loop.last %}, {% endif %}{%endfor%}) const {
    {% for i in range(N) %}
        if(i{{i}}<i{{i}}_start|i{{i}}>i{{i}}_stop){            
            printf("%s: index i{{i}}=%ld out of range (%ld, %ld)\\n", array_name.c_str(), i{{i}}, i{{i}}_start, i{{i}}_stop);  
            exit(EXIT_FAILURE);
        }
    {%endfor%}
    }
#endif
    
    inline const T& operator()({% for i in range(N) %}int i{{i}}{% if not loop.last %}, {% endif %}{%endfor%}) const {
#ifdef MULTIARRAY_INDICES_CHECK
        check_indices({% for i in range(N) %}i{{i}}{% if not loop.last %}, {% endif %}{%endfor%});
#endif        
        return data[{% for i in range(N) %} (i{{i}}-i{{i}}_start){% if not loop.last %}*s[{{i}}]+{% endif %}{%endfor%}];
    }
    inline T& operator()({% for i in range(N) %}int i{{i}}{% if not loop.last %}, {% endif %}{%endfor%}) {
#ifdef MULTIARRAY_INDICES_CHECK    
        check_indices({% for i in range(N) %}i{{i}}{% if not loop.last %}, {% endif %}{%endfor%});
#endif        
        return data[{% for i in range(N) %} (i{{i}}-i{{i}}_start){% if not loop.last %}*s[{{i}}]+{% endif %}{%endfor%}]; 
    }
};
"""
general_class_template = Template(cpp_general_class_code_template)


# dimension
def generate_template_class(N):
    return class_template.render(N=N)


def generate_general_template_class(N):
    return general_class_template.render(N=N)


header_file_template = Template("""
//automatically generate source code

#ifndef ACE_MULTIARRAY_H
#define ACE_MULTIARRAY_H
#include <string.h>
#include "contigous_array_nd.h"

using namespace std;

{% for class_template in class_templates %}
{{class_template}}
{% endfor %}


#endif //ACE_MULTIARRAY_H
""")

class_templates = [generate_template_class(N) for N in range(1, DIMMAX)]
genral_class_templates = [generate_general_template_class(N) for N in range(1, DIMMAX)]
print(header_file_template.render(class_templates=class_templates + genral_class_templates))
