"""Deals with the socket communication between the PIMD and driver code.

Copyright (C) 2013, Joshua More and Michele Ceriotti

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <http.//www.gnu.org/licenses/>.


Deals with creating the socket, transmitting and receiving data, accepting and
removing different driver routines and the parallelization of the force
calculation.

Classes:
   Status: Simple class to keep track of the status, uses bitwise or to give
      combinations of different status options.
   DriverSocket: Class to deal with communication between a client and
      the driver code.
   InterfaceSocket: Host server class. Deals with distribution of all the jobs
      between the different client servers.

Functions:
   Message: Sends a header string through the socket.

Exceptions:
   Disconnected: Raised if client has been disconnected.
   InvalidStatus: Raised if client has the wrong status. Shouldn't have to be
      used if the structure of the program is correct.
"""

__all__ = ['InterfaceSocket']

import numpy as np
import sys, os
import socket, select, threading, signal, string, time
from ipi.utils.depend import depstrip
from ipi.utils.messages import verbosity, warning, info
from ipi.utils.softexit import softexit


HDRLEN = 12
UPDATEFREQ = 10
TIMEOUT = 5.0
SERVERTIMEOUT = 2.0*TIMEOUT
NTIMEOUT = 10

def Message(mystr):
   """Returns a header of standard length HDRLEN."""

   return string.ljust(string.upper(mystr), HDRLEN)


class Disconnected(Exception):
   """Disconnected: Raised if client has been disconnected."""

   pass

class InvalidSize(Exception):
   """Disconnected: Raised if client returns forces with inconsistent number of atoms."""

   pass

class InvalidStatus(Exception):
   """InvalidStatus: Raised if client has the wrong status.

   Shouldn't have to be used if the structure of the program is correct.
   """

   pass

class Status:
   """Simple class used to keep track of the status of the client.

   Uses bitwise or to give combinations of different status options.
   i.e. Status.Up | Status.Ready would be understood to mean that the client
   was connected and ready to receive the position and cell data.

   Attributes:
      Disconnected: Flag for if the client has disconnected.
      Up: Flag for if the client is running.
      Ready: Flag for if the client has ready to receive position and cell data.
      NeedsInit: Flag for if the client is ready to receive forcefield
         parameters.
      HasData: Flag for if the client is ready to send force data.
      Busy: Flag for if the client is busy.
      Timeout: Flag for if the connection has timed out.
   """

   Disconnected = 0
   Up = 1
   Ready = 2
   NeedsInit = 4
   HasData = 8
   Busy = 16
   Timeout = 32


class DriverSocket(socket.socket):
   """Deals with communication between the client and driver code.

   Deals with sending and receiving the data from the driver code. Keeps track
   of the status of the driver. Initialises the driver forcefield, sends the
   position and cell data, and receives the force data.

   Attributes:
      _buf: A string buffer to hold the reply from the driver.
      status: Keeps track of the status of the driver.
      lastreq: The ID of the last request processed by the client.
      locked: Flag to mark if the client has been working consistently on one image.
   """

   def __init__(self, socket):
      """Initialises DriverSocket.

      Args:
         socket: A socket through which the communication should be done.
      """

      super(DriverSocket,self).__init__(_sock=socket)
      self._buf = np.zeros(0,np.byte)
      self.peername = self.getpeername()
      self.status = Status.Up
      self.waitstatus = False
      self.lastreq = None
      self.locked = False
      
   def shutdown(self, how=socket.SHUT_RDWR):
      
      self.sendall(Message("exit"))
      self.status = Status.Disconnected
      super(DriverSocket,self).shutdown(how)
   
   def poll(self):
      """Waits for driver status."""

      self.status = Status.Disconnected  # sets disconnected as failsafe status, in case _getstatus fails and exceptions are ignored upstream
      self.status = self._getstatus()

   def _getstatus(self):
      """Gets driver status.

      Returns:
         An integer labelling the status via bitwise or of the relevant members
         of Status.
      """
            
      if not self.waitstatus:
         try:
            readable, writable, errored = select.select([], [self], [])
            if self in writable:
               self.sendall(Message("status"))
               self.waitstatus = True
         except:
            return Status.Disconnected

      try:
         reply = self.recv(HDRLEN)
         self.waitstatus = False # got status reply         
      except socket.timeout:
         warning(" @SOCKET:   Timeout in status recv!", verbosity.debug )
         return Status.Up | Status.Busy | Status.Timeout
      except:
         return Status.Disconnected
      
      if not len(reply) == HDRLEN:
         return Status.Disconnected
      elif reply == Message("ready"):
         return Status.Up | Status.Ready
      elif reply == Message("needinit"):
         return Status.Up | Status.NeedsInit
      elif reply == Message("havedata"):
         return Status.Up | Status.HasData
      else:
         warning(" @SOCKET:    Unrecognized reply: " + str(reply), verbosity.low )
         return Status.Up

   def recvall(self, dest):
      """Gets the potential energy, force and virial from the driver.

      Args:
         dest: Object to be read into.

      Raises:
         Disconnected: Raised if client is disconnected.

      Returns:
         The data read from the socket to be read into dest.
      """

      blen = dest.itemsize*dest.size
      if (blen > len(self._buf)):
         self._buf.resize(blen)
      bpos = 0
      ntimeout = 0

      while bpos < blen:
         timeout = False

#   pre-2.5 version.
         try:
            bpart = ""            
            bpart = self.recv(blen - bpos)
            if len(bpart) == 0: raise socket.timeout  # There is a problem if this returns no data
            self._buf[bpos:bpos + len(bpart)] = np.fromstring(bpart, np.byte)
         except socket.timeout:
            warning(" @SOCKET:   Timeout in status recvall, trying again!", verbosity.low)
            timeout = True
            ntimeout += 1
            if ntimeout > NTIMEOUT:
               warning(" @SOCKET:  Couldn't receive within %5d attempts. Time to give up!" % (NTIMEOUT), verbosity.low)
               raise Disconnected()
            pass
         if (not timeout and bpart == 0):
            raise Disconnected()
         bpos += len(bpart)

#   post-2.5 version: slightly more compact for modern python versions
#         try:
#            bpart = 1
#            bpart = self.recv_into(self._buf[bpos:], blen-bpos)
#         except socket.timeout:
#            print " @SOCKET:   Timeout in status recvall, trying again!"
#            timeout = True
#            pass
#         if (not timeout and bpart == 0):
#            raise Disconnected()
#         bpos += bpart
#TODO this Disconnected() exception currently just causes the program to hang.
#This should do something more graceful

      if np.isscalar(dest):
         return np.fromstring(self._buf[0:blen], dest.dtype)[0]
      else:
         return np.fromstring(self._buf[0:blen], dest.dtype).reshape(dest.shape)

   def initialize(self, rid, pars):
      """Sends the initialization string to the driver.

      Args:
         rid: The index of the request, i.e. the replica that
            the force calculation is for.
         pars: The parameter string to be sent to the driver.

      Raises:
         InvalidStatus: Raised if the status is not NeedsInit.
      """

      if self.status & Status.NeedsInit:
         try:
            self.sendall(Message("init"))
            self.sendall(np.int32(rid))
            self.sendall(np.int32(len(pars)))
            self.sendall(pars)
         except:
            self.poll()
            return
      else:
         raise InvalidStatus("Status in init was " + self.status)

   def sendpos(self, pos, cell):
      """Sends the position and cell data to the driver.

      Args:
         pos: An array containing the atom positions.
         cell: A cell object giving the system box.

      Raises:
         InvalidStatus: Raised if the status is not Ready.
      """

      if (self.status & Status.Ready):
         try:
            self.sendall(Message("posdata"))
            self.sendall(cell.h, 9*8)
            self.sendall(cell.ih, 9*8)
            self.sendall(np.int32(len(pos)/3))
            self.sendall(pos, len(pos)*8)
         except:
            self.poll()
            return
      else:
         raise InvalidStatus("Status in sendpos was " + self.status)

   def getforce(self):
      """Gets the potential energy, force and virial from the driver.

      Raises:
         InvalidStatus: Raised if the status is not HasData.
         Disconnected: Raised if the driver has disconnected.

      Returns:
         A list of the form [potential, force, virial, extra].
      """

      if (self.status & Status.HasData):
         self.sendall(Message("getforce"));
         reply = ""         
         while True:            
            try:
               reply = self.recv(HDRLEN)
            except socket.timeout:
               warning(" @SOCKET:   Timeout in getforce, trying again!", verbosity.low)
               continue
            if reply == Message("forceready"):
               break
            else:
               warning(" @SOCKET:   Unexpected getforce reply: %s" % (reply), verbosity.low)
            if reply == "":
               raise Disconnected()
      else:
         raise InvalidStatus("Status in getforce was " + self.status)
      
      mu = np.float64()
      mu = self.recvall(mu)

      mlen = np.int32()
      mlen = self.recvall(mlen)
      mf = np.zeros(3*mlen,np.float64)
      mf = self.recvall(mf)

      mvir = np.zeros((3,3),np.float64)
      mvir = self.recvall(mvir)
      
      #! Machinery to return a string as an "extra" field. Comment if you are using a old patched driver that does not return anything!
      mlen = np.int32()
      mlen = self.recvall(mlen)
      if mlen > 0 :
         mxtra = np.zeros(mlen,np.character)
         mxtra = self.recvall(mxtra)
         mxtra = "".join(mxtra)
      else:
         mxtra = ""

      #!TODO must set up a machinery to intercept the "extra" return field
      return [mu, mf, mvir, mxtra]


class InterfaceSocket(object):
   """Host server class.

   Deals with distribution of all the jobs between the different client servers
   and both initially and as clients either finish or are disconnected.
   Deals with cleaning up after all calculations are done. Also deals with the
   threading mechanism, and cleaning up if the interface is killed.

   Attributes:
      address: A string giving the name of the host network.
      port: An integer giving the port the socket will be using.
      slots: An integer giving the maximum allowed backlog of queued clients.
      mode: A string giving the type of socket used.
      latency: A float giving the number of seconds the interface will wait
         before updating the client list.
      timeout: A float giving a timeout limit for considering a calculation dead
         and dropping the connection.
      dopbc: A boolean which decides whether or not to fold the bead positions
         back into the unit cell before passing them to the client code.
      server: The socket used for data transmition.
      clients: A list of the driver clients connected to the server.
      requests: A list of all the jobs required in the current PIMD step.
      jobs: A list of all the jobs currently running.
      _poll_thread: The thread the poll loop is running on.
      _prev_kill: Holds the signals to be sent to clean up the main thread
         when a kill signal is sent.
      _poll_true: A boolean giving whether the thread is alive.
      _poll_iter: An integer used to decide whether or not to check for
         client connections. It is used as a counter, once it becomes higher
         than the pre-defined number of steps between checks the socket will
         update the list of clients and then be reset to zero.
   """

   def __init__(self, address="localhost", port=31415, slots=4, mode="unix", latency=1e-3, timeout=1.0, dopbc=True):
      """Initialises interface.

      Args:
         address: An optional string giving the name of the host server.
            Defaults to 'localhost'.
         port: An optional integer giving the port number. Defaults to 31415.
         slots: An optional integer giving the maximum allowed backlog of
            queueing clients. Defaults to 4.
         mode: An optional string giving the type of socket. Defaults to 'unix'.
         latency: An optional float giving the time in seconds the socket will
            wait before updating the client list. Defaults to 1e-3.
         timeout: Length of time waiting for data from a client before we assume
            the connection is dead and disconnect the client.
         dopbc: A boolean which decides whether or not to fold the bead positions
            back into the unit cell before passing them to the client code.

      Raises:
         NameError: Raised if mode is not 'unix' or 'inet'.
      """

      self.address = address
      self.port = port
      self.slots = slots
      self.mode = mode
      self.latency = latency
      self.timeout = timeout
      self.dopbc = dopbc
      self._poll_thread = None
      self._prev_kill = {}
      self._poll_true = False
      self._poll_iter = 0

   def open(self):
      """Creates a new socket.

      Used so that we can create a interface object without having to also
      create the associated socket object.
      """

      if self.mode == "unix":
         self.server = socket.socket(socket.AF_UNIX, socket.SOCK_STREAM)
         try:
            self.server.bind("/tmp/ipi_" + self.address)
            info("Created unix socket with address " + self.address, verbosity.medium)
         except:
            raise ValueError("Error opening unix socket. Check if a file " + ("/tmp/ipi_" + self.address) + " exists, and remove it if unused.")

      elif self.mode == "inet":
         self.server = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
         self.server.bind((self.address,self.port))
         info("Created inet socket with address " + self.address + " and port number " + str(self.port), verbosity.medium)
      else:
         raise NameError("InterfaceSocket mode " + self.mode + " is not implemented (should be unix/inet)")

      self.server.listen(self.slots)
      self.server.settimeout(SERVERTIMEOUT)
      self.clients = []
      self.requests = []
      self.jobs = []

   def close(self):
      """Closes down the socket."""

      info(" @SOCKET: Shutting down the driver interface.", verbosity.low )
      
      for c in self.clients[:]:
         if (c.status & Status.Up):
            c.shutdown(socket.SHUT_RDWR)
            
      self.server.shutdown(socket.SHUT_RDWR)
      self.server.close()
      if self.mode == "unix":
         os.unlink("/tmp/ipi_" + self.address)

   def queue(self, atoms, cell, pars=None, reqid=0):
      """Adds a request.

      Note that the pars dictionary need to be sent as a string of a
      standard format so that the initialization of the driver can be done.

      Args:
         atoms: An Atoms object giving the atom positions.
         cell: A Cell object giving the system box.
         pars: An optional dictionary giving the parameters to be sent to the
            driver for initialization. Defaults to {}.
         reqid: An optional integer that identifies requests of the same type,
            e.g. the bead index

      Returns:
         A list giving the status of the request of the form {'atoms': Atoms
         object giving the atom positions, 'cell': Cell object giving the
         system box, 'pars': parameter string, 'result': holds the result as a
         list once the computation is done, 'status': a string labelling the
         status, 'id': the id of the request, usually the bead number, 'start':
         the starting time for the calculation, used to check for timeouts.}.
      """

      par_str = " "

      if not pars is None:
         for k,v in pars.items():
            par_str += k + " : " + str(v) + " , "
      else:
         par_str = " "

      # APPLY PBC -- this is useful for codes such as LAMMPS that don't do full PBC when computing distances
      pbcpos = depstrip(atoms.q).copy()
      if self.dopbc:
         cell.array_pbc(pbcpos)

      newreq = {"pos": pbcpos, "cell": cell, "pars": par_str,
                "result": None, "status": "Queued", "id": reqid,
                "start": -1 }

      self.requests.append(newreq)
      return newreq

   def release(self, request):
      """Empties the list of requests once finished.

      Args:
         request: A list of requests that are done.
      """

      if request in self.requests:
         self.requests.remove(request)

   def pool_update(self):
      """Deals with keeping the pool of client drivers up-to-date during a
      force calculation step.

      Deals with maintaining the client list. Clients that have
      disconnected are removed and their jobs removed from the list of
      running jobs and new clients are connected to the server.
      """

      for c in self.clients[:]:
         if not (c.status & Status.Up):
            try:
               warning(" @SOCKET:   Client " + str(c.peername) +" died or got unresponsive(C). Removing from the list.", verbosity.low)
               c.shutdown(socket.SHUT_RDWR)
               c.close()
            except:
               pass
            c.status = Status.Disconnected
            self.clients.remove(c)
            for [k,j] in self.jobs[:]:
               if j is c:
                  self.jobs = [ w for w in self.jobs if not ( w[0] is k and w[1] is j ) ] # removes pair in a robust way
                  #self.jobs.remove([k,j])
                  k["status"] = "Queued"
                  k["start"] = -1

      keepsearch = True
      while keepsearch:
         readable, writable, errored = select.select([self.server], [], [], 0.0)
         if self.server in readable:
            client, address = self.server.accept()
            client.settimeout(TIMEOUT)
            driver = DriverSocket(client)
            info(" @SOCKET:   Client asked for connection from "+ str( address ) +". Now hand-shaking.", verbosity.low)
            driver.poll()
            if (driver.status | Status.Up):
               self.clients.append(driver)
               info(" @SOCKET:   Handshaking was successful. Added to the client list.", verbosity.low)
            else:
               warning(" @SOCKET:   Handshaking failed. Dropping connection.", verbosity.low)
               client.shutdown(socket.SHUT_RDWR)
               client.close()
         else:
            keepsearch = False

   def pool_distribute(self):
      """Deals with keeping the list of jobs up-to-date during a force
      calculation step.

      Deals with maintaining the jobs list. Gets data from drivers that have
      finished their calculation and removes that job from the list of running
      jobs, adds jobs to free clients and initialises the forcefields of new
      clients.
      """

      for c in self.clients:
         if c.status == Status.Disconnected : # client disconnected. force a pool_update
            self._poll_iter = UPDATEFREQ
            return
         if not c.status & ( Status.Ready | Status.NeedsInit ):
            c.poll()

      for [r,c] in self.jobs[:]:
         if c.status & Status.HasData:
            try:
               r["result"] = c.getforce()
               if len(r["result"][1]) != len(r["pos"]):
                  raise InvalidSize
            except Disconnected:
               c.status = Status.Disconnected
               continue
            except InvalidSize:
              warning(" @SOCKET:   Client returned an inconsistent number of forces. Will mark as disconnected and try to carry on.", verbosity.low)
              c.status = Status.Disconnected
              continue
            except:
              warning(" @SOCKET:   Client got in a awkward state during getforce. Will mark as disconnected and try to carry on.", verbosity.low)
              c.status = Status.Disconnected
              continue
            c.poll()
            while c.status & Status.Busy: # waits, but check if we got stuck.
               if self.timeout > 0 and r["start"] > 0 and time.time() - r["start"] > self.timeout:
                  warning(" @SOCKET:  Timeout! HASDATA for bead " + str(r["id"]) + " has been running for " + str(time.time() - r["start"]) + " sec.", verbosity.low)
                  warning(" @SOCKET:   Client " + str(c.peername) + " died or got unresponsive(A). Disconnecting.", verbosity.low)
                  try:
                     c.shutdown(socket.SHUT_RDWR)
                  except:
                     pass
                  c.close()
                  c.status = Status.Disconnected
                  continue
               c.poll()
            if not (c.status & Status.Up):
               warning(" @SOCKET:   Client died a horrible death while getting forces. Will try to cleanup.", verbosity.low)
               continue
            r["status"] = "Done"
            c.lastreq = r["id"] # saves the ID of the request that the client has just processed
            self.jobs = [ w for w in self.jobs if not ( w[0] is r and w[1] is c ) ] # removes pair in a robust way

         if self.timeout > 0 and c.status != Status.Disconnected and r["start"] > 0 and time.time() - r["start"] > self.timeout:
            warning(" @SOCKET:  Timeout! Request for bead " + str( r["id"]) + " has been running for " + str(time.time() - r["start"]) + " sec.", verbosity.low)
            warning(" @SOCKET:   Client " + str(c.peername) + " died or got unresponsive(B). Disconnecting.",verbosity.low)
            try:
               c.shutdown(socket.SHUT_RDWR)
            except socket.error:
               e = sys.exc_info()
               warning(" @SOCKET:  could not shut down cleanly the socket. %s: %s in file '%s' on line %d" % (e[0].__name__, e[1], os.path.basename(e[2].tb_frame.f_code.co_filename), e[2].tb_lineno), verbosity.low )
            c.close()
            c.poll()
            c.status = Status.Disconnected

      freec = self.clients[:]
      for [r2, c] in self.jobs:
         freec.remove(c)

      pendr = self.requests[:]
      pendr = [ r for r in self.requests if r["status"] == "Queued" ]

      for fc in freec[:]:
         matched = False
         # first, makes sure that the client is REALLY free
         if not (fc.status & Status.Up):
            self.clients.remove(fc)   # if fc is in freec it can't be associated with a job (we just checked for that above)
            continue
         if fc.status & Status.HasData:
            continue
         if not (fc.status & (Status.Ready | Status.NeedsInit | Status.Busy) ):
            warning(" @SOCKET: Client " + str(fc.peername) + " is in an unexpected status " + str(fc.status) + " at (1). Will try to keep calm and carry on.", verbosity.low)
            continue
         for match_ids in ( "match", "none", "free", "any" ):
            for r in pendr[:]:
               if match_ids == "match" and not fc.lastreq is r["id"]:
                  continue
               elif match_ids == "none" and not fc.lastreq is None:
                  continue
               elif match_ids == "free" and fc.locked:
                  continue

               info(" @SOCKET: Assigning [%5s] request id %4s to client with last-id %4s (% 3d/% 3d : %s)" % (match_ids,  str(r["id"]),  str(fc.lastreq), self.clients.index(fc), len(self.clients), str(fc.peername) ), verbosity.high )

               while fc.status & Status.Busy:
                  fc.poll()
               if fc.status & Status.NeedsInit:
                  fc.initialize(r["id"], r["pars"])
                  fc.poll()
                  while fc.status & Status.Busy: # waits for initialization to finish. hopefully this is fast
                     fc.poll()
               if fc.status & Status.Ready:
                  fc.sendpos(r["pos"], r["cell"])
                  r["status"] = "Running"
                  r["start"] = time.time() # sets start time for the request
                  fc.poll()
                  self.jobs.append([r,fc])
                  fc.locked =  (fc.lastreq is r["id"])
                  matched = True
                  # removes r from the list of pending jobs
                  pendr = [nr for nr in pendr if (not nr is r)]
                  break
               else:
                  warning(" @SOCKET: Client " + str(fc.peername) + " is in an unexpected status " + str(fc.status) + " at (2). Will try to keep calm and carry on.", verbosity.low)
            if matched:
               break # doesn't do a second (or third) round if it managed
                     # to assign the job

   def _kill_handler(self, signal, frame):
      """Deals with handling a kill call gracefully.

      Prevents any of the threads becoming zombies, by intercepting a
      kill signal using the standard python function signal.signal() and
      then closing the socket and the spawned threads before closing the main
      thread. Called when signals SIG_INT and SIG_TERM are received.

      Args:
         signal: An integer giving the signal number of the signal received
            from the socket.
         frame: Current stack frame.
      """

      warning(" @SOCKET:   Kill signal. Trying to make a clean exit.", verbosity.low)
      self.end_thread()

      softexit.trigger(" @SOCKET: Kill signal received")

      try:
         self.__del__()
      except:
         pass
      if signal in self._prev_kill:
         self._prev_kill[signal](signal, frame)

   def _poll_loop(self):
      """The main thread loop.

      Runs until either the program finishes or a kill call is sent. Updates
      the pool of clients every UPDATEFREQ loops and loops every latency
      seconds until _poll_true becomes false.
      """

      info(" @SOCKET: Starting the polling thread main loop.", verbosity.low)
      self._poll_iter = UPDATEFREQ
      while self._poll_true:
         time.sleep(self.latency)
         # makes sure to remove the last dead client as soon as possible -- and to get clients if we are dry
         if self._poll_iter >= UPDATEFREQ or len(self.clients)==0 or (len(self.clients) > 0 and not(self.clients[0].status & Status.Up)):
            self.pool_update()
            self._poll_iter = 0
         self._poll_iter += 1
         self.pool_distribute()

         if os.path.exists("EXIT"): # softexit
            info(" @SOCKET: Soft exit request from file EXIT. Flushing job queue.", verbosity.low)
            # releases all pending requests
            for r in self.requests:
               r["status"] = "Exit"
            for c in self.clients:
               try:
                  c.shutdown(socket.SHUT_RDWR)
                  c.close()
               except:
                  pass
            # flush it all down the drain
            self.clients = []
            self.jobs = []
      self._poll_thread = None

   def started(self):
      """Returns a boolean specifying whether the thread has started yet."""

      return (not self._poll_thread is None)

   def start_thread(self):
      """Spawns a new thread.

      Splits the main program into two threads, one that runs the polling loop
      which updates the client list, and one which gets the data. Also sets up
      the machinery to deal with a kill call, in the case of a Ctrl-C or
      similar signal the signal is intercepted by the _kill_handler function,
      which cleans up the spawned thread before closing the main thread.

      Raises:
         NameError: Raised if the polling thread already exists.
      """

      self.open()
      if not self._poll_thread is None:
         raise NameError("Polling thread already started")
      self._poll_thread = threading.Thread(target=self._poll_loop, name="poll_" + self.address)
      self._poll_thread.daemon = True
      self._prev_kill[signal.SIGINT] = signal.signal(signal.SIGINT, self._kill_handler)
      self._prev_kill[signal.SIGTERM] = signal.signal(signal.SIGTERM, self._kill_handler)
      self._poll_true = True
      self._poll_thread.start()

   def end_thread(self):
      """Closes the spawned thread.

      Deals with cleaning up the spawned thread cleanly. First sets
      _poll_true to false to indicate that the poll_loop should be exited, then
      closes the spawned thread and removes it.
      """

      self._poll_true = False
      if not self._poll_thread is None:
         self._poll_thread.join()
      self._poll_thread = None
      self.close()
