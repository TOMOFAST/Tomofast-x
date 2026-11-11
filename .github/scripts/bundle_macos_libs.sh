#!/bin/bash
#
# Copies every non-system dynamic library that a macOS executable needs into a
# "lib" folder next to it, and rewrites the load paths so that the executable
# finds them there. The result runs on a Mac without Homebrew or gfortran.
#
# Usage: bundle_macos_libs.sh <bundle_dir> <executable_name> <fallback_lib_dir>...
#   bundle_dir        folder that holds the executable (lib/ is created inside it)
#   executable_name   file name of the executable inside bundle_dir
#   fallback_lib_dir  folders searched when a library is referenced via @rpath
#
set -euo pipefail

BUNDLE_DIR="$1"
EXE="$BUNDLE_DIR/$2"
shift 2
SEARCH_DIRS=("$@")

LIB_DIR="$BUNDLE_DIR/lib"
mkdir -p "$LIB_DIR"

# Libraries that ship with macOS and must not be bundled.
is_system_lib() {
  case "$1" in
    /usr/lib/*|/System/*) return 0 ;;
    *) return 1 ;;
  esac
}

# Prints the dependencies of a Mach-O file (skips the first line and the file's own id).
list_deps() {
  local file="$1"
  local own_id
  own_id=$(otool -D "$file" | tail -n +2 || true)
  otool -L "$file" | tail -n +2 | awk '{print $1}' | while read -r dep; do
    [ "$dep" = "$own_id" ] && continue
    echo "$dep"
  done
}

# Finds the real file behind a dependency path (resolves @rpath via the search dirs).
resolve_dep() {
  local dep="$1"
  local name
  name=$(basename "$dep")
  case "$dep" in
    @rpath/*|@loader_path/*|@executable_path/*)
      for dir in "${SEARCH_DIRS[@]}"; do
        if [ -f "$dir/$name" ]; then
          echo "$dir/$name"
          return 0
        fi
      done
      echo "ERROR: cannot find $dep in: ${SEARCH_DIRS[*]}" >&2
      return 1
      ;;
    *)
      echo "$dep"
      ;;
  esac
}

# Copies the dependencies of one file into lib/ and points the file at the copies.
process_file() {
  local file="$1"
  local dep real name
  for dep in $(list_deps "$file"); do
    is_system_lib "$dep" && continue
    name=$(basename "$dep")
    if [ ! -f "$LIB_DIR/$name" ]; then
      real=$(resolve_dep "$dep")
      echo "  bundling $name (from $real)"
      cp -L "$real" "$LIB_DIR/$name"
      chmod u+w "$LIB_DIR/$name"
      install_name_tool -id "@rpath/$name" "$LIB_DIR/$name"
      # Recurse into the new library's own dependencies.
      process_file "$LIB_DIR/$name"
    fi
    install_name_tool -change "$dep" "@rpath/$name" "$file"
  done
}

echo "Bundling libraries for $EXE"
process_file "$EXE"

# Remove any rpaths pointing at the build machine, then add the bundle-relative ones.
strip_rpaths() {
  local file="$1"
  otool -l "$file" | awk '/cmd LC_RPATH/{getline; getline; print $2}' | while read -r rp; do
    install_name_tool -delete_rpath "$rp" "$file" 2>/dev/null || true
  done
}

strip_rpaths "$EXE"
install_name_tool -add_rpath "@executable_path/lib" "$EXE"
for lib in "$LIB_DIR"/*.dylib; do
  strip_rpaths "$lib"
  install_name_tool -add_rpath "@loader_path" "$lib"
done

# install_name_tool invalidates signatures; Apple Silicon refuses to run unsigned code,
# so re-sign everything with an ad-hoc signature.
for f in "$LIB_DIR"/*.dylib "$EXE"; do
  codesign --force --sign - "$f"
done

# Fail if anything still points outside the bundle or the system.
echo "Checking the bundle for external references..."
bad=0
for f in "$EXE" "$LIB_DIR"/*.dylib; do
  for dep in $(list_deps "$f"); do
    is_system_lib "$dep" && continue
    case "$dep" in
      @rpath/*)
        [ -f "$LIB_DIR/$(basename "$dep")" ] || { echo "MISSING: $f -> $dep"; bad=1; } ;;
      *)
        echo "EXTERNAL: $f -> $dep"; bad=1 ;;
    esac
  done
done
[ "$bad" -eq 0 ] || exit 1

echo "Bundle contents:"
ls -l "$LIB_DIR"
otool -L "$EXE"
