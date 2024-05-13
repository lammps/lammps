/*
 *_________________________________________________________________________*
 *      POEMS: PARALLELIZABLE OPEN SOURCE EFFICIENT MULTIBODY SOFTWARE     *
 *      DESCRIPTION: SEE READ-ME                                           *
 *      FILE NAME: poemslist.h                                             *
 *      AUTHORS: See Author List                                           *
 *      GRANTS: See Grants List                                            *
 *      COPYRIGHT: (C) 2005 by Authors as listed in Author's List          *
 *      LICENSE: Please see License Agreement                              *
 *      DOWNLOAD: Free at www.rpi.edu/~anderk5                             *
 *      ADMINISTRATOR: Prof. Kurt Anderson                                 *
 *                     Computational Dynamics Lab                          *
 *                     Rensselaer Polytechnic Institute                    *
 *                     110 8th St. Troy NY 12180                           *
 *      CONTACT:        anderk5@rpi.edu                                    *
 *_________________________________________________________________________*/

#ifndef _POEMSLIST_H_
#define _POEMSLIST_H_

#include <iostream>
#include <cstdlib>

using namespace std;

template<class T> class ListElement{
public:
  ListElement<T>* prev;
  ListElement<T>* next;
  T* value;

  ListElement();
  ListElement(T* v);
  ~ListElement();
};

template<class S> class List {
  int numelements;
  ListElement<S>* head;
  ListElement<S>* tail;
public:
  List();
  ~List();
  int GetNumElements();
  ListElement<S>* GetHeadElement();
  ListElement<S>* GetTailElement();
  void Remove(ListElement<S>* ele);
  ListElement<S>* Append(S* v);
  ListElement<S>* Prepend(S* v);
  S** CreateArray();
  S* operator()(int id);
  void Append(List<S> * listToAppend);
  void DeleteValues();
  void RemoveElementAndDeleteValue(ListElement<S>* ele);
  void PrintList();
};

//
// ListElement
//

template<class T> ListElement<T>::ListElement(){
  next = prev = 0;
  value = 0;
}

template<class T> ListElement<T>::ListElement(T* v){
  next = prev = nullptr;
  value = v;
}

template<class T> ListElement<T>::~ListElement(){
}

//
// List
//

template<class S> List<S>::List(){
  head = tail = nullptr;
  numelements = 0;
}

template<class S> List<S>::~List(){
  // delete all list elements

  while(numelements)
    Remove(tail);
}

template<class S> void List<S>::Append(List<S> * listToAppend)
{
        tail->next = listToAppend->head;
        listToAppend->head->prev = tail;
        tail = listToAppend->tail;
}

template<class S> int List<S>::GetNumElements(){
  return numelements;
}

template<class S> ListElement<S>* List<S>::GetHeadElement(){
  return head;
}

template<class S> ListElement<S>* List<S>::GetTailElement(){
        return tail;
}


template<class S> void List<S>::Remove(ListElement<S>* ele){
  if(!ele){
    cerr << "ERROR: ListElement to be removed not defined" << endl;
    exit(0);
  }
  if(!numelements){
    cerr << "ERROR: List is empty" << endl;
    exit(0);
  }

  if(ele != head)
    ele->prev->next = ele->next;
  else
    head = ele->next;

  if(ele != tail)
    ele->next->prev = ele->prev;
  else
    tail = ele->prev;

  numelements--;
  if(ele)
  delete ele;
}

template<class S> ListElement<S>* List<S>::Append(S* v){
  if(!v){
    cerr << "ERROR: cannot add null Body to list" << endl;
    exit(0);
  }

  numelements++;
  ListElement<S>* ele = new ListElement<S>(v);

  if(numelements==1)
    head = tail = ele;
  else{
          /*
    tail->next = ele;
    ele->prev = tail;
         tail = ele;*/

          ele->prev = tail;
          tail = ele;
          ele->prev->next = ele;

  }
  return ele;
}

template<class S> ListElement<S>* List<S>::Prepend(S* v){
  if(!v){
    cerr << "ERROR: cannot add null Body to list" << endl;
    exit(0);
  }

  numelements++;
  ListElement<S>* ele = new ListElement<S>(v);

  if(numelements==1)
    head = tail = ele;
  else{
          ele->next = head;
          head = ele;
          ele->next->prev = ele;
  }
  return ele;
}

template<class S> S** List<S>::CreateArray(){
  S** array = new S* [numelements];

  ListElement<S>* ele = head;
  for(int i=0;ele != nullptr;i++){
    array[i] = ele->value;
    ele = ele->next;
  }
  return array;
}

template<class S> S* List<S>::operator()(int id){
  if(id >= numelements){
    cerr << "ERROR: subscript out of bounds" << endl;
    exit(0);
  }

  ListElement<S>* ele = head;
  for(int i=0;i<id;i++){
    ele = ele->next;
  }

  return ele->value;
}

template<class S> void List<S>::DeleteValues(){
  while(numelements)
    RemoveElementAndDeleteValue(tail);
}

template<class S> void List<S>::RemoveElementAndDeleteValue(ListElement<S>* ele){
  S* v = ele->value;
  Remove(ele);
  delete v;
}

template<class S> void List<S>::PrintList(){
        cout<<"Printing List "<<endl;
        ListElement<S>* ele = head;
        cout<<*(ele->value)<<" ";
        ele = ele->next;
        for(int k =2; k<numelements; k++){
                cout<<*(ele->value)<<" ";
                ele = ele->next;
        }
        cout<<*(ele->value)<<endl;
}


#endif
