# MCP tools

[Deep Agents Code](/oss/python/deepagents/code/overview)

# MCP tools
Copy page

Load additional tools from MCP (Model Context Protocol) serversCopy page
> 

## Documentation Index

Fetch the complete documentation index at: [https://docs.langchain.com/llms.txt](https://docs.langchain.com/llms.txt)

Use this file to discover all available pages before exploring further.[MCP (Model Context Protocol)](https://modelcontextprotocol.io/) lets you extend Deep Agents Code with tools from external servers—file systems, APIs, databases, and more—without modifying the agent itself. Deep Agents Code connects to MCP servers at startup, discovers their tools, and makes them available to the agent alongside the built-in tools.
Add MCP servers by adding a `.mcp.json` config file to your project for project-level scope, or at user-level to apply to all projects.

## [​](#quickstart)Quickstart

This quickstart adds the [LangChain documentation MCP server](https://docs.langchain.com/mcp) to every Deep Agents Code session on your machine. Swap in any other MCP server’s URL or stdio command in the same shape.
[](#)

Create the config fileIf not already present, create the `.mcp.json` file at user-level to make the server available to every project on the machine or at project-level.
{' ' * (self.list_depth - 1)}- User
{' ' * (self.list_depth - 1)}- Project
{' ' * (self.list_depth - 1)}- Project (hidden)
```
mkdir -p ~/.deepagents
touch ~/.deepagents/.mcp.json

```
Servers in this file (`~/.deepagents/.mcp.json`) are available in every project on this machine.
```
touch .mcp.json

```
Servers in this file (`<project>/.mcp.json`) are available to this project.
```
mkdir -p .deepagents
touch .deepagents/.mcp.json

```
Servers in this file (`<project>/.deepagents/.mcp.json`) are available to this project but kept out of the repo root.See [Discovery locations](#discovery-locations) for full precedence rules.[](#)

Add the MCP server~/.deepagents/.mcp.json
```
{
 "mcpServers": {
 "docs-langchain": {
 "type": "http",
 "url": "https://docs.langchain.com/mcp"
 }
 }
}

```
To add more servers, add more entries to `mcpServers`. See [Configuration format](#configuration-format) for OAuth, stdio, SSE, and HTTP server fields, environment variables, and headers.[](#)

Launch Deep Agents Code
```
dcode

```
On startup, Deep Agents Code auto-discovers the config, connects to each server, discovers its tools, and prints a confirmation:
```
✓ Loaded 3 MCP tools

```
Run `/mcp` in an interactive session to see per-server status, transport, and the loaded tool list. The agent can now use those tools for the duration of the session—stdio servers are kept alive between tool calls.

## [​](#auto-discovery)Auto-discovery

Deep Agents Code automatically searches for `.mcp.json` files in standard locations. No flags are needed—just place a config file and it gets picked up.

### [​](#discovery-locations)Discovery locations

Configs are checked in this order (lowest to highest precedence):

| 
 | Priority | Location | Scope
 | 1 (lowest) | `~/.deepagents/.mcp.json` | User-level—applies to all projects
 | 2 | `<project>/.deepagents/.mcp.json` | Project-level—`.deepagents` subdirectory
 | 3 (highest) | `<project>/.mcp.json` | Project-level—root (Claude Code compatible)
The project root is the nearest parent directory containing a `.git` folder, falling back to the current working directory.
When multiple config files exist, their `mcpServers` entries are merged. If the same server name appears in more than one file, the higher-precedence config wins. This lets a project-level config override a user-level entry (for example, pinning a different version of the same server) without disturbing your other projects.

### [​](#flags)Flags

| 
 | Flag | Behavior
 | `--mcp-config PATH` | Add an explicit config as the highest-precedence source (merged on top of auto-discovered configs)
 | `--no-mcp` | Disable MCP entirely—no servers are loaded
`--mcp-config` and `--no-mcp` are mutually exclusive.

### [​](#claude-code-compatibility)Claude Code compatibility

If you already have a `.mcp.json` at your project root for Claude Code, Deep Agents Code picks it up automatically—no extra setup needed.

## [​](#configuration-format)Configuration format

Each key under `mcpServers` is a server name. The server’s fields determine how Deep Agents Code connects to it.

### [​](#stdio-servers-default)stdio servers (default)

stdio servers are spawned as child processes. Deep Agents Code communicates with them over stdin/stdout.
mcp-config.json
```
{
 "mcpServers": {
 "filesystem": {
 "command": "npx",
 "args": ["-y", "@modelcontextprotocol/server-filesystem", "/tmp"],
 "env": {}
 },
 "github": {
 "command": "npx",
 "args": ["-y", "@modelcontextprotocol/server-github"],
 "env": { "GITHUB_TOKEN": "your-token" }
 }
 }
}

```

### [​](#sse-and-http-servers)SSE and HTTP servers

For remote MCP servers, set `type` to `"sse"` or `"http"` and provide a `url`:
mcp-config.json
```
{
 "mcpServers": {
 "remote-api": {
 "type": "sse",
 "url": "https://api.example.com/mcp",
 "headers": { "Authorization": "Bearer your-token" }
 }
 }
}

```

### [​](#field-reference)Field reference

stdio (default)**Required:** `command`. **Optional:** `args`, `env`, plus the shared [tool-filter fields](#tool-filtering).[​](#param-command)commandstringrequiredThe executable to run.[​](#param-args)argsstring[]Arguments passed to the command.[​](#param-env)envobjectEnvironment variables set for the subprocess. Use this to pass API keys and other credentials without exposing them in shell history.

sse**Required:** `type: "sse"`, `url`. **Optional:** `headers`, `auth`, plus the shared [tool-filter fields](#tool-filtering).[​](#param-type)type"sse"requiredTransport type. Use `"sse"` for Server-Sent Events.[​](#param-url)urlstringrequiredThe server endpoint URL.[​](#param-headers)headersobjectHTTP headers sent with every request. Commonly used for authentication. Values support `${VAR}` references to parent-shell environment variables (resolved when the server activates).[​](#param-auth)auth"oauth"Set to `"oauth"` to drive an OAuth login flow with `dcode mcp login` instead of supplying an `Authorization` header. Cannot be combined with an `Authorization` header. See [OAuth login](#oauth-login).

http**Required:** `type: "http"`, `url`. **Optional:** `headers`, `auth`, plus the shared [tool-filter fields](#tool-filtering).[​](#param-type-1)type"http"requiredTransport type. Use `"http"` for streamable HTTP. `streamable_http` and `streamable-http` are accepted as aliases.[​](#param-url-1)urlstringrequiredThe server endpoint URL.[​](#param-headers-1)headersobjectHTTP headers sent with every request. Commonly used for authentication. Values support `${VAR}` references to parent-shell environment variables (resolved when the server activates).[​](#param-auth-1)auth"oauth"Set to `"oauth"` to drive an OAuth login flow with `dcode mcp login` instead of supplying an `Authorization` header. Cannot be combined with an `Authorization` header. See [OAuth login](#oauth-login).
The `type` field can also be written as `transport` for compatibility with other MCP clients.
Server names must match `[A-Za-z0-9_-]+`. Names are used as on-disk basenames for OAuth token files, so path separators and other shell metacharacters are rejected at config load.

### [​](#header-environment-variables)Header environment variables

Header values support `${VAR}` substitution from the parent shell, resolved at server activation rather than at config load. One unset variable only fails the server that needs it; the rest still come up.
.mcp.json
```
{
 "mcpServers": {
 "internal-api": {
 "type": "http",
 "url": "https://api.example.com/mcp",
 "headers": { "Authorization": "Bearer ${INTERNAL_API_TOKEN}" }
 }
 }
}

```

## [​](#multiple-servers)Multiple servers

You can configure as many servers as you need. Tools from all servers are merged and available to the agent:
mcp-config.json
```
{
 "mcpServers": {
 "filesystem": {
 "command": "npx",
 "args": ["-y", "@modelcontextprotocol/server-filesystem", "/home/user/projects"]
 },
 "github": {
 "command": "npx",
 "args": ["-y", "@modelcontextprotocol/server-github"],
 "env": { "GITHUB_TOKEN": "ghp_..." }
 },
 "database": {
 "type": "sse",
 "url": "https://db-mcp.internal:8080/mcp",
 "headers": { "Authorization": "Bearer ..." }
 }
 }
}

```

## [​](#tool-filtering)Tool filtering

Each server may narrow the tools it exposes to the agent with one of two optional fields:

{' ' * (self.list_depth - 1)}- `allowedTools`: keep only the listed tools; drop everything else.

{' ' * (self.list_depth - 1)}- `disabledTools`: drop the listed tools; keep everything else.

Filtering applies to stdio, HTTP, and SSE servers alike. Both of the following are rejected at config load:

{' ' * (self.list_depth - 1)}- Setting `allowedTools` and `disabledTools` on the same server.

{' ' * (self.list_depth - 1)}- Setting either field to an empty list (would silently strip every tool, or be a no-op). Omit the field instead.

.mcp.json
```
{
 "mcpServers": {
 "filesystem": {
 "command": "npx",
 "args": ["-y", "@modelcontextprotocol/server-filesystem", "/tmp"],
 "allowedTools": ["read_file", "list_directory"]
 },
 "github": {
 "command": "npx",
 "args": ["-y", "@modelcontextprotocol/server-github"],
 "disabledTools": ["delete_repository", "delete_*_branch"]
 }
 }
}

```

### [​](#match-rules)Match rules

Each entry is a literal tool name or an [`fnmatch`](https://docs.python.org/3/library/fnmatch.html)-style glob (any entry containing `*`, `?`, or `[` is treated as a pattern). Entries are matched against both the bare MCP tool name and the server-prefixed form (`{server}_{tool}`), so either form works:

```
{
 "allowedTools": ["read_file", "fs_list_*"]
}

```

Entries that match no loaded tool are logged as a warning, not an error — the underlying MCP server can evolve its tool list across versions without breaking your config.
[​](#param-allowed-tools)allowedToolsstring[]Tool names or `fnmatch` glob patterns to keep. All other tools from this server are dropped. Mutually exclusive with `disabledTools`.
[​](#param-disabled-tools)disabledToolsstring[]Tool names or `fnmatch` glob patterns to drop. All other tools from this server are kept. Mutually exclusive with `allowedTools`.

## [​](#oauth-login)OAuth login

For remote MCP servers that require OAuth (Slack, GitHub, Notion, Linear, and other hosted MCP endpoints), set `"auth": "oauth"` on the server entry and run the login subcommand once. Tokens are persisted to disk and refreshed automatically.

### [​](#configure-the-server)Configure the server

.mcp.json
```
{
 "mcpServers": {
 "linear": {
 "type": "http",
 "url": "https://mcp.linear.app/mcp",
 "auth": "oauth"
 }
 }
}

```

`auth: "oauth"` is mutually exclusive with an `Authorization` header on the same entry, and cannot be set on a stdio server.

### [​](#run-the-login-flow)Run the login flow

```
dcode mcp login linear

```

What happens depends on the server’s host:

{' ' * (self.list_depth - 1)}- **Spec-compliant servers** (the default): Deep Agents Code performs Dynamic Client Registration, opens an Authorization Code + PKCE flow in your browser, and asks you to paste the redirected URL back into the terminal.

{' ' * (self.list_depth - 1)}- **Slack** (`slack.com`, `*.slack.com`): same paste-back flow, but with Slack’s public client preseeded. You’re prompted for an optional team ID (e.g., `T01234567`) so the app installs into the right workspace.

{' ' * (self.list_depth - 1)}- **GitHub** (`api.githubcopilot.com`): RFC 8628 Device Authorization Grant. Deep Agents Code prints a verification URL and a user code; you enter the code in your browser and Deep Agents Code polls for completion.

By default, `dcode mcp login` reads the same auto-discovered configs Deep Agents Code uses at runtime (subject to project-level trust gating). Pass `--config <path>` to use a specific file:

```
dcode mcp login linear --config ./mcp-config.json

```

Project-level configs that haven’t been trusted (see [Project-level trust](#project-level-trust)) are skipped during `mcp login` to prevent attacker-controlled `headers` entries from exfiltrating local secrets through `${VAR}` interpolation. Run `dcode` in the project once to approve the config, or pass `--config <path>` explicitly.

### [​](#token-storage)Token storage

Tokens are written to:

```
~/.deepagents/.state/mcp-tokens/<server>-<sha256-16(url)>.json

```

The `<sha256-16(url)>` segment is the first 16 hex characters of the SHA-256 of the server URL. The directory is locked to mode `0700` and each token file is mode `0600`. Files include the OAuth access token, refresh token, and the dynamically registered client info, all in a schema-versioned payload that’s written atomically (write-to-temp + `rename`).
Hashing the URL into the filename means the same server name pointing at different URLs (for example, dev vs. prod) gets independent token files and can’t trample each other.

### [​](#re-authentication)Re-authentication

When refresh fails at runtime (the refresh token expired or was revoked), Deep Agents Code marks the server as `unauthenticated` instead of crashing the agent. The welcome banner shows the count of unauthenticated servers, and `/mcp` reports the reason per server. Re-run `dcode mcp login <server>` to refresh credentials — your conversation continues without restarting.

## [​](#server-status)Server status

Each configured server lands in one of three states after startup:

| 
 | Status | Meaning
 | `ok` | Connected; tools are loaded and available to the agent
 | `unauthenticated` | OAuth login required or refresh failed — run `dcode mcp login <server>`
 | `error` | Pre-flight, discovery, or transport setup failed; an error message is attached
A single failing server no longer aborts startup. The agent runs with whichever servers came up cleanly, and the welcome banner surfaces counts of unauthenticated and errored servers next to the tool count. Open `/mcp` in an interactive session to see per-server status, transport, tool list, and the failure reason for non-`ok` entries. The viewer live-updates as servers connect and supports `tab`/`shift+tab` navigation.

## [​](#project-level-trust)Project-level trust

Project-level configs can contain stdio servers that execute local commands and remote servers whose `headers` may interpolate `${VAR}` from your environment. To prevent untrusted repositories from running arbitrary code or exfiltrating local secrets on CLI startup, Deep Agents Code enforces a **default-deny** policy for project-level entries.

### [​](#how-it-works)How it works

{' ' * (self.list_depth - 1)}- **Interactive mode:** Deep Agents Code prompts for approval before activating project servers, showing each stdio command and remote URL. Approval is persisted using a SHA-256 content fingerprint—if the config changes, you are prompted again.

{' ' * (self.list_depth - 1)}- **Non-interactive mode (`-n`):** Project servers are silently skipped unless `--trust-project-mcp` is passed.

{' ' * (self.list_depth - 1)}- **Trust covers stdio and remote entries alike** — remote servers can SSRF into localhost or cloud-metadata endpoints during the pre-flight probe and exfiltrate `${VAR}` values via headers, so they’re gated the same way as stdio.

{' ' * (self.list_depth - 1)}- **User-level configs** (`~/.deepagents/.mcp.json`) are always trusted—the same trust model as `config.toml` and `hooks.json`.

{' ' * (self.list_depth - 1)}- **`dcode mcp login`** also honors project trust: an untrusted project-level config is skipped during login discovery so an attacker-controlled remote entry cannot pull secrets into the OAuth handshake.

### [​](#flags-2)Flags

| 
 | Flag | Behavior
 | `--trust-project-mcp` | Trust all project-level stdio servers without prompting (for CI and automation)

```
# Skip the approval prompt
dcode --trust-project-mcp

# Non-interactive: explicitly trust project servers
dcode -n "run tests" --trust-project-mcp

```

### [​](#trust-store)Trust store

Trust decisions are stored in `~/.deepagents/config.toml`:

```
[mcp_trust.projects]
"/Users/you/myproject" = "sha256:abc123..."

```

Each key is an absolute project root path. The value is a SHA-256 digest of the concatenated project-level config contents. To revoke trust, delete the entry or modify the project’s `.mcp.json` (which invalidates the fingerprint automatically).
A trusted stdio MCP server has the same permissions as your user account. Only approve servers from repositories you trust. Review the commands shown in the approval prompt before accepting.

## [​](#system-prompt-awareness)System prompt awareness

Connected MCP servers and their tools are automatically listed in the agent’s system prompt, grouped by server name and transport type. This helps the model reason about tool provenance and failure domains without requiring manual context.

## [​](#troubleshooting)Troubleshooting

Server fails to start (stdio)Verify the command works outside Deep Agents Code:
```
npx -y @modelcontextprotocol/server-filesystem /tmp

```
Common causes: the package isn’t installed, `npx` isn’t on `PATH`, or required environment variables are missing.

Connection refused (SSE/HTTP)Check that the remote server is running and the URL is correct. If the server requires authentication, make sure `headers` includes the correct credentials.

Tools not appearingDeep Agents Code prints the number of tools loaded at startup (e.g., `✓ Loaded 3 MCP tools`). If you see `0`, the server started successfully but didn’t advertise any tools—check the server’s own logs or documentation.

Server shows `unauthenticated` in /mcpEither you haven’t run `dcode mcp login <server>` yet, or the persisted refresh token expired or was revoked server-side. Run the login command again — your session keeps running and the server will re-attach once tokens are refreshed.

`Invalid MCP config at ...`A pre-flight validation rejected `--mcp-config` (or an auto-discovered `.mcp.json`). Common causes: an unsupported server name (must match `[A-Za-z0-9_-]+`), `auth: oauth` on a stdio server, both `command` and `url` set on the same entry, or a header value that isn’t a string. Fix the highlighted reason and relaunch — Deep Agents Code no longer dumps a multi-page subprocess trace for config errors.

`${VAR}` header references failHeader interpolation runs at activation time, so an unset variable only fails the server that needs it. Export the variable in the parent shell or add it to `~/.deepagents/.env`. To debug, set `DEEPAGENTS_CODE_DEBUG=1` and inspect the per-session log path printed to stderr on shutdown.

## [​](#further-reading)Further reading

{' ' * (self.list_depth - 1)}- [LangChain MCP guide](/oss/python/langchain/mcp): protocol details, building custom servers, and using `langchain-mcp-adapters` programmatically

{' ' * (self.list_depth - 1)}- [MCP specification](https://modelcontextprotocol.io/): the official protocol spec and server registry

---

[Connect these docs](/use-these-docs) to Claude, VSCode, and more via MCP for real-time answers.[Edit this page on GitHub](https://github.com/langchain-ai/docs/edit/main/src/oss/deepagents/code/mcp-tools.mdx) or [file an issue](https://github.com/langchain-ai/docs/issues/new/choose).

Was this page helpful?YesNo[ConfigurationPrevious](/oss/python/deepagents/code/configuration)[Data locationsNext](/oss/python/deepagents/code/data-locations)⌘I[Docs by LangChain home page
![light logo](https://mintcdn.com/langchain-5e9cc07a/nQm-sjd_MByLhgeW/images/brand/langchain-docs-dark-blue.png?fit=max&auto=format&n=nQm-sjd_MByLhgeW&q=85&s=5babf1a1962208fd7eed942fa2432ecb)
![dark logo](https://mintcdn.com/langchain-5e9cc07a/nQm-sjd_MByLhgeW/images/brand/langchain-docs-light-blue.png?fit=max&auto=format&n=nQm-sjd_MByLhgeW&q=85&s=0bcd2a1f2599ed228bcedf0f535b45b1)](/)[github](https://github.com/langchain-ai)[x](https://x.com/LangChain)[linkedin](https://www.linkedin.com/company/langchain)[youtube](https://www.youtube.com/@LangChain)

Resources[Forum](https://forum.langchain.com/)[Changelog](https://changelog.langchain.com/)[LangChain Academy](https://academy.langchain.com/)[Contact Sales](https://www.langchain.com/contact-sales)

Company[Home](https://langchain.com/)[Trust Center](https://trust.langchain.com/)[Careers](https://langchain.com/careers)[Blog](https://blog.langchain.com/)[github](https://github.com/langchain-ai)[x](https://x.com/LangChain)[linkedin](https://www.linkedin.com/company/langchain)[youtube](https://www.youtube.com/@LangChain)
