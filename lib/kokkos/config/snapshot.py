#! /usr/bin/env python

"""
Snapshot a project into another project and perform the necessary repo actions
to provide a commit message that can be used to trace back to the exact point
in the source repository.
"""

#todo:
#  Support svn
#  Allow renaming of the source dir in the destination path
#  Check if a new snapshot is necessary?
#

import sys

#check the version number so that there is a good error message when argparse is not available.
#This checks for exactly 2.7 which is bad, but it is a python 2 script and argparse was introduced
#in 2.7 which is also the last version of python 2. If this script is updated for python 3 this
#will need to change, but for now it is not safe to allow 3.x to run this.
if sys.version_info[:2] != (2, 7):
  print "Error snapshot requires python 2.7 detected version is %d.%d." % (sys.version_info[0], sys.version_info[1])
  sys.exit(1)

import subprocess, argparse, re, doctest, os, datetime, traceback

def parse_cmdline(description):
  parser = argparse.ArgumentParser(usage="snapshot.py [options] source destination", description=description)

  parser.add_argument("-n", "--no-commit", action="store_false", dest="create_commit", default=True,
                      help="Do not perform a commit or create a commit message.")
  parser.add_argument("-v", "--verbose", action="store_true", dest="verbose_mode", default=False,
                      help="Enable verbose mode.")
  parser.add_argument("-d", "--debug", action="store_true", dest="debug_mode", default=False,
                      help="Enable debugging output.")
  parser.add_argument("--no-validate-repo", action="store_true", dest="no_validate_repo", default=False,
                      help="Reduce the validation that the source and destination repos are clean to a warning.")
  parser.add_argument("--source-repo", choices=["git","none"], default="",
                      help="Type of repository of the source, use none to skip all repository operations.")
  parser.add_argument("--dest-repo", choices=["git","none"], default="",
                      help="Type of repository of the destination, use none to skip all repository operations.")
  parser.add_argument("--small", action="store_true", dest="small_mode",
                      help="Don't include tests and other extra files when copying.")

  parser.add_argument("source",      help="Source project to snapshot from.")
  parser.add_argument("destination", help="Destination to snapshot too.")

  options = parser.parse_args()
  options = validate_options(options)
  return options
#end parseCmdline

def validate_options(options):
  apparent_source_repo_type="none"
  apparent_dest_repo_type="none"

  #prevent user from accidentally giving us a path that rsync will treat differently than expected.
  options.source      = options.source.rstrip(os.sep)
  options.destination = options.destination.rstrip(os.sep)

  options.source      = os.path.abspath(options.source)
  options.destination = os.path.abspath(options.destination)

  if os.path.exists(options.source):
    apparent_source_repo_type, source_root = determine_repo_type(options.source)
  else:
    raise RuntimeError("Could not find source directory of %s." % options.source)
  options.source_root = source_root

  if not os.path.exists(options.destination):
    print "Could not find destination directory of %s so it will be created." % options.destination
    os.makedirs(options.destination)

  apparent_dest_repo_type, dest_root = determine_repo_type(options.destination)
  options.dest_root = dest_root

  #error on svn repo types for now
  if apparent_source_repo_type == "svn" or apparent_dest_repo_type == "svn":
    raise RuntimeError("SVN repositories are not supported at this time.")

  if options.source_repo == "":
    #source repo type is not specified to just using the apparent type.
    options.source_repo = apparent_source_repo_type
  else:
    if options.source_repo != "none" and options.source_repo != apparent_source_repo_type:
      raise RuntimeError("Specified source repository type of %s conflicts with determined type of %s" % \
        (options.source_repo, apparent_source_repo_type))

  if options.dest_repo == "":
    #destination repo type is not specified to just using the apparent type.
    options.dest_repo = apparent_dest_repo_type
  else:
    if options.dest_repo != "none" and options.dest_repo != apparent_dest_repo_type:
      raise RuntimeError("Specified destination repository type of %s conflicts with determined type of %s" % \
        (options.dest_repo, apparent_dest_repo_type))

  return options
#end validate_options

def run_cmd(cmd, options, working_dir="."):
  cmd_str = " ".join(cmd)
  if options.verbose_mode:
    print "Running command '%s' in dir %s." % (cmd_str, working_dir)

  proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=working_dir)
  proc_stdout, proc_stderr = proc.communicate()
  ret_val = proc.wait()

  if options.debug_mode:
    print "==== %s stdout start ====" % cmd_str
    print proc_stdout
    print "==== %s stdout end ====" % cmd_str
    print "==== %s stderr ====" % cmd_str
    print proc_stderr
    print "==== %s stderr ====" % cmd_str

  if ret_val != 0:
    raise RuntimeError("Command '%s' failed with error code %d. Error message:%s%s%sstdout:%s" % \
      (cmd_str, ret_val, os.linesep, proc_stderr, os.linesep, proc_stdout))

  return proc_stdout, proc_stderr
#end run_cmd

def determine_repo_type(location):
  apparent_repo_type = "none"

  while location != "":
    if os.path.exists(os.path.join(location, ".git")):
      apparent_repo_type = "git"
      break
    elif os.path.exists(os.path.join(location, ".svn")):
      apparent_repo_type = "svn"
      break
    else:
      location = location[:location.rfind(os.sep)]

  return apparent_repo_type, location
#end determine_repo_type

def rsync(source, dest, options):
  rsync_cmd = ["rsync", "-ar", "--delete"]
  if options.debug_mode:
    rsync_cmd.append("-v")

  if options.small_mode or options.source_repo == "git":
    rsync_cmd.append("--delete-excluded")

  if options.small_mode:
    rsync_cmd.append("--include=config/master_history.txt")
    rsync_cmd.append("--include=cmake/tpls")
    rsync_cmd.append("--exclude=benchmarks/")
    rsync_cmd.append("--exclude=config/*")
    rsync_cmd.append("--exclude=doc/")
    rsync_cmd.append("--exclude=example/")
    rsync_cmd.append("--exclude=tpls/")
    rsync_cmd.append("--exclude=HOW_TO_SNAPSHOT")
    rsync_cmd.append("--exclude=unit_test")
    rsync_cmd.append("--exclude=unit_tests")
    rsync_cmd.append("--exclude=perf_test")
    rsync_cmd.append("--exclude=performance_tests")

  if options.source_repo == "git":
    rsync_cmd.append("--exclude=.git*")

  rsync_cmd.append(options.source)
  rsync_cmd.append(options.destination)
  run_cmd(rsync_cmd, options)
#end rsync

def create_commit_message(commit_id, commit_log, project_name, project_location):
  eol = os.linesep
  message = "Snapshot of %s from commit %s" % (project_name, commit_id)
  message += eol * 2
  message += "From repository at %s" % project_location
  message += eol * 2
  message += "At commit:" + eol
  message += commit_log
  return message
#end create_commit_message

def find_git_commit_information(options):
  r"""
  >>> class fake_options:
  ...   source="."
  ...   verbose_mode=False
  ...   debug_mode=False
  >>> myoptions = fake_options()
  >>> find_git_commit_information(myoptions)[2:]
  ('sems', 'software.sandia.gov:/git/sems')
  """
  git_log_cmd = ["git", "log", "-1"]

  output, error = run_cmd(git_log_cmd, options, options.source)

  commit_match = re.match("commit ([0-9a-fA-F]+)", output)
  commit_id = commit_match.group(1)
  commit_log = output

  git_remote_cmd = ["git", "remote", "-v"]
  output, error = run_cmd(git_remote_cmd, options, options.source)

  remote_match = re.search("origin\s([^ ]*/([^ ]+))", output, re.MULTILINE)
  if not remote_match:
    raise RuntimeError("Could not find origin of repo at %s. Consider using none for source repo type." % (options.source))

  source_location = remote_match.group(1)
  source_name     = remote_match.group(2).strip()

  if source_name[-1] == "/":
    source_name = source_name[:-1]

  return commit_id, commit_log, source_name, source_location
#end find_git_commit_information

def do_git_commit(message, options):
  if options.verbose_mode:
    print "Commiting to destination repository."

  git_add_cmd = ["git", "add", "-A"]
  run_cmd(git_add_cmd, options, options.destination)

  git_commit_cmd = ["git", "commit", "-m%s" % message]
  run_cmd(git_commit_cmd, options, options.destination)

  git_log_cmd = ["git", "log", "--format=%h", "-1"]
  commit_sha1, error = run_cmd(git_log_cmd, options, options.destination)

  print "Commit %s was made to %s." % (commit_sha1.strip(), options.dest_root)
#end do_git_commit

def verify_git_repo_clean(location, options):
  git_status_cmd = ["git", "status", "--porcelain"]
  output, error = run_cmd(git_status_cmd, options, location)

  if output != "":
    if options.no_validate_repo == False:
      raise RuntimeError("%s is not clean.%sPlease commit or stash all changes before running snapshot."
        % (location, os.linesep))
    else:
      print "WARNING: %s is not clean. Proceeding anyway." % location
      print "WARNING:   This could lead to differences in the source and destination."
      print "WARNING:   It could also lead to extra files being included in the snapshot commit."
#end verify_git_repo_clean

def main(options):
  if options.verbose_mode:
    print "Snapshotting %s to %s." % (options.source, options.destination)

  if options.source_repo == "git":
    verify_git_repo_clean(options.source, options)
    commit_id, commit_log, repo_name, repo_location = find_git_commit_information(options)
  elif options.source_repo == "none":
    commit_id     = "N/A"
    commit_log    = "Unknown commit from %s snapshotted at: %s" % (options.source, datetime.datetime.now())
    repo_name     = options.source
    repo_location = options.source

  commit_message = create_commit_message(commit_id, commit_log, repo_name, repo_location) + os.linesep*2

  if options.dest_repo == "git":
    verify_git_repo_clean(options.destination, options)

  rsync(options.source, options.destination, options)

  if options.dest_repo == "git":
    do_git_commit(commit_message, options)
  elif options.dest_repo == "none":
    file_name = "snapshot_message.txt"
    message_file = open(file_name, "w")
    message_file.write(commit_message)
    message_file.close()
    cwd = os.getcwd()
    print "No commit done by request. Please use file at:"
    print "%s%sif you wish to commit this to a repo later." % (cwd+"/"+file_name, os.linesep)
#end main

if (__name__ == "__main__"):
  if ("--test" in sys.argv):
    doctest.testmod()
    sys.exit(0)

  try:
    options = parse_cmdline(__doc__)
    main(options)
  except RuntimeError, e:
    print "Error occured:", e
    if "--debug" in sys.argv:
      traceback.print_exc()
    sys.exit(1)
  else:
    sys.exit(0)
